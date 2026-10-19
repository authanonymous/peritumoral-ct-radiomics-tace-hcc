"""
Peritumoral Radiomics + Habitat Imaging — 3D Slicer Script
============================================================
Paste into 3D Slicer Python Interactor (Ctrl+3), or run via
exec(open("<path to this file>").read()).

Notes on implementation choices:
  - Local std (_fast_local_std) uses the identity Var(x) = E[x^2] - E[x]^2
    via nd.uniform_filter rather than nd.generic_filter(arr, np.std, ...),
    which is substantially faster on full-volume CT arrays.
  - All SimpleITK/ITK operations run on the main thread; SimpleITK is not
    thread-safe, so only pure NumPy/sklearn work (k-means, radiomics
    extraction) runs on the worker thread.
  - distance_transform_edt is computed once; all shell bands are sliced
    from the single resulting distance array.
  - resampledPixelSpacing=None avoids double-resampling inside PyRadiomics.
  - peripheral_enhancement_ratio uses the absolute HU denominator (not the
    liver-mean-normalised value).
"""

import slicer, sitkUtils, os, datetime, traceback
import threading
import queue as _queue
import numpy as np
import SimpleITK as sitk
import scipy.ndimage as nd
import qt

try:
    from radiomics import featureextractor
except ImportError:
    slicer.util.pip_install("pyradiomics")
    from radiomics import featureextractor

try:
    from sklearn.cluster import KMeans
except ImportError:
    slicer.util.pip_install("scikit-learn")
    from sklearn.cluster import KMeans

import pandas as pd


# ======================================================================
#  SECTION 1 — SCENE HELPERS
# ======================================================================

def get_volume_names():
    names = []
    col = slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode")
    col.InitTraversal()
    n = col.GetNextItemAsObject()
    while n:
        names.append(n.GetName())
        n = col.GetNextItemAsObject()
    return names

def get_mask_names():
    entries = []
    col = slicer.mrmlScene.GetNodesByClass("vtkMRMLLabelMapVolumeNode")
    col.InitTraversal()
    n = col.GetNextItemAsObject()
    while n:
        entries.append(n.GetName())
        n = col.GetNextItemAsObject()
    col = slicer.mrmlScene.GetNodesByClass("vtkMRMLSegmentationNode")
    col.InitTraversal()
    n = col.GetNextItemAsObject()
    while n:
        seg = n.GetSegmentation()
        for i in range(seg.GetNumberOfSegments()):
            sid   = seg.GetNthSegmentID(i)
            sname = seg.GetSegment(sid).GetName()
            entries.append("%s > %s" % (n.GetName(), sname))
        n = col.GetNextItemAsObject()
    return entries

def pull_mask_as_sitk(entry, ct_node):
    if " > " in entry:
        node_name, seg_name = entry.split(" > ", 1)
        seg_node = slicer.mrmlScene.GetFirstNodeByName(node_name)
        if seg_node is None:
            raise ValueError("Segmentation node not found: " + node_name)
        seg    = seg_node.GetSegmentation()
        seg_id = None
        for i in range(seg.GetNumberOfSegments()):
            sid = seg.GetNthSegmentID(i)
            if seg.GetSegment(sid).GetName() == seg_name:
                seg_id = sid
                break
        if seg_id is None:
            raise ValueError("Segment '%s' not found in '%s'" % (seg_name, node_name))
        tmp = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLLabelMapVolumeNode", "__tmp_mask__")
        slicer.modules.segmentations.logic().ExportSegmentsToLabelmapNode(
            seg_node, [seg_id], tmp, ct_node)
        img = sitkUtils.PullVolumeFromSlicer(tmp)
        slicer.mrmlScene.RemoveNode(tmp)
    else:
        node = slicer.mrmlScene.GetFirstNodeByName(entry)
        if node is None:
            raise ValueError("Mask node not found: " + entry)
        img = sitkUtils.PullVolumeFromSlicer(node)
    return sitk.Cast(img, sitk.sitkUInt8)


# ======================================================================
#  SECTION 2 — IMAGE PROCESSING
# ======================================================================

SHELLS = [
    ("peri_0-3mm",  0.0,  3.0,  (0.20, 0.78, 0.20)),
    ("peri_3-5mm",  3.0,  5.0,  (1.00, 0.60, 0.00)),
    ("peri_5-10mm", 5.0, 10.0,  (0.85, 0.12, 0.12)),
]

HABITAT_COLOURS = {
    "habitat_1_low":   (0.20, 0.40, 0.85),
    "habitat_2_mid":   (0.95, 0.77, 0.06),
    "habitat_3_high":  (0.85, 0.12, 0.12),
    "habitat_4_extra": (0.55, 0.10, 0.80),
}

RIM_COLOURS = {
    "Tumor (reference)": (0.22, 0.42, 0.88),
    "peri_0-3mm":        (0.20, 0.78, 0.20),
    "peri_3-5mm":        (1.00, 0.60, 0.00),
    "peri_5-10mm":       (0.85, 0.12, 0.12),
}

def resample_to_1mm(img):
    sp     = img.GetSpacing()
    sz     = img.GetSize()
    new_sz = [int(round(sz[i] * sp[i])) for i in range(3)]
    f = sitk.ResampleImageFilter()
    f.SetOutputSpacing([1.0, 1.0, 1.0])
    f.SetSize(new_sz)
    f.SetOutputDirection(img.GetDirection())
    f.SetOutputOrigin(img.GetOrigin())
    f.SetInterpolator(sitk.sitkLinear)
    return f.Execute(img)

def resample_mask(mask, ref):
    return sitk.BinaryThreshold(
        sitk.Resample(mask, ref, sitk.Transform(),
                      sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8),
        1, 1000, 1, 0)

def make_shells_fast(t_mask, l_mask, shells_on, logfn=print):
    """Distance transform is computed once; all shell bands are sliced from it."""
    t_arr = sitk.GetArrayFromImage(t_mask).astype(bool)
    l_arr = sitk.GetArrayFromImage(l_mask).astype(bool)

    logfn("  Computing distance transform (once for all shells)...")
    spacing = t_mask.GetSpacing()[::-1]   # sitk x,y,z → numpy z,y,x
    dist    = nd.distance_transform_edt(~t_arr, sampling=spacing)

    shell_masks = {}
    for name, lo, hi, _ in SHELLS:
        if name not in shells_on:
            continue
        band = (dist > lo) & (dist <= hi) & l_arr
        nvx  = int(np.sum(band))
        if nvx < 50:
            logfn("  %s : SKIPPED (%d vox)" % (name, nvx))
            continue
        img = sitk.GetImageFromArray(band.astype(np.uint8))
        img.CopyInformation(t_mask)
        shell_masks[name] = img
        logfn("  %s : %d vox  OK" % (name, nvx))

    return shell_masks

def make_extractor(bin_width, shape=True, fast=False):
    """fast=True restricts extraction to first-order features only (used for small habitat masks)."""
    e = featureextractor.RadiomicsFeatureExtractor()
    e.disableAllFeatures()
    classes = ['firstorder'] if fast else ['firstorder', 'glcm', 'glszm', 'glrlm', 'gldm']
    for c in classes:
        e.enableFeatureClassByName(c)
    if shape:
        e.enableFeatureClassByName('shape')
    e.settings.update({
        'binWidth':              bin_width,
        'resampledPixelSpacing': None,
        'randomSeed':            42,
        'correctMask':           True,
    })
    return e


# ======================================================================
#  SECTION 3 — HABITAT IMAGING
# ======================================================================

def _fast_local_std(arr):
    """3x3x3 local standard deviation via nd.uniform_filter (Var(x) = E[x^2] - E[x]^2)."""
    arr_f   = arr.astype(np.float32)
    ex      = nd.uniform_filter(arr_f,    size=3)
    ex2     = nd.uniform_filter(arr_f**2, size=3)
    var     = np.maximum(ex2 - ex**2, 0.0)   # clamp floating-point negatives
    return np.sqrt(var)

def compute_habitats(ct_norm, t_mask, n_habitats=3,
                     use_texture=True, logfn=print):
    """
    Partition tumour voxels into n_habitats using k-means.
    Features per voxel:
      - normalised HU intensity
      - local 3x3x3 std
      - gradient magnitude
    """
    ct_arr = sitk.GetArrayFromImage(ct_norm)
    t_arr  = sitk.GetArrayFromImage(t_mask).astype(bool)
    indices   = np.argwhere(t_arr)
    intensity = ct_arr[t_arr].reshape(-1, 1)

    feature_vecs = [intensity]

    if use_texture:
        logfn("  Computing local std (fast)...")
        local_std = _fast_local_std(ct_arr)
        feature_vecs.append(local_std[t_arr].reshape(-1, 1))

        logfn("  Computing gradient magnitude...")
        arr_f32  = ct_arr.astype(np.float32)
        gz, gy, gx = np.gradient(arr_f32)
        grad_mag = np.sqrt(gz**2 + gy**2 + gx**2, dtype=np.float32)
        feature_vecs.append(grad_mag[t_arr].reshape(-1, 1))
        del gz, gy, gx, grad_mag, arr_f32   # free memory immediately

    X = np.hstack(feature_vecs).astype(np.float32)

    # Normalise each feature to [0, 1]
    X_min = X.min(axis=0)
    X_rng = X.max(axis=0) - X_min
    X_rng[X_rng < 1e-8] = 1.0
    X_norm = (X - X_min) / X_rng

    logfn("  Running k-means (k=%d, %d voxels)..." % (n_habitats, len(X_norm)))
    km     = KMeans(n_clusters=n_habitats, random_state=42, n_init=10)
    labels = km.fit_predict(X_norm)

    # Sort clusters by mean intensity (low → high) for stable naming
    cluster_means = [ct_arr[t_arr][labels == c].mean() for c in range(n_habitats)]
    order         = np.argsort(cluster_means)
    remap         = {old: new for new, old in enumerate(order)}
    labels_sorted = np.array([remap[l] for l in labels])

    label_vol = np.zeros_like(ct_arr, dtype=np.int32)
    for idx_3d, lab in zip(indices, labels_sorted):
        z, y, x = idx_3d
        label_vol[z, y, x] = lab + 1

    habitat_masks = {}
    cluster_stats = {}
    name_map      = ["low", "mid", "high", "extra"]

    for h in range(n_habitats):
        hname = "habitat_%d_%s" % (h + 1, name_map[h] if h < len(name_map) else str(h + 1))
        hmask = (label_vol == h + 1).astype(np.uint8)
        nvox  = int(np.sum(hmask))
        hmean = float(ct_arr[hmask.astype(bool)].mean()) if nvox > 0 else 0.0

        sitk_h = sitk.GetImageFromArray(hmask)
        sitk_h.CopyInformation(t_mask)
        habitat_masks[hname] = sitk_h
        cluster_stats[hname] = {
            "voxels":  nvox,
            "mean_hu": hmean,
            "vol_mm3": nvox * float(np.prod(t_mask.GetSpacing())),
        }
        logfn("  %s: %d vox | mean HU %.1f" % (hname, nvox, hmean))

    return habitat_masks, cluster_stats


def compute_habitat_metrics(cluster_stats, n_habitats, ct_norm, t_mask):
    total_vox = sum(s["voxels"] for s in cluster_stats.values())
    if total_vox == 0:
        return {}

    fractions = [s["voxels"] / total_vox for s in cluster_stats.values()]
    metrics   = {}

    fracs_nz = [f for f in fractions if f > 0]
    metrics["habitat_entropy"]   = float(-np.sum([f * np.log2(f) for f in fracs_nz]))
    metrics["habitat_dominance"] = float(max(fractions))

    for h, (hname, stats) in enumerate(cluster_stats.items()):
        metrics["habitat_%d_vol_fraction" % (h + 1)] = stats["voxels"] / total_vox
        metrics["habitat_%d_mean_hu"      % (h + 1)] = stats["mean_hu"]

    viable_vox = sum(s["voxels"] for k, s in cluster_stats.items()
                     if not k.endswith("_low"))
    metrics["viable_fraction"]     = float(viable_vox / total_vox)
    means = [s["mean_hu"] for s in cluster_stats.values()]
    metrics["inter_habitat_contrast"] = float(max(means) - min(means))

    low_key = [k for k in cluster_stats if k.endswith("_low")]
    if low_key:
        metrics["necrotic_fraction"] = cluster_stats[low_key[0]]["voxels"] / total_vox

    return metrics


# ======================================================================
#  SECTION 4 — VISUALISATION & SAVE
# ======================================================================

def show_rims_in_scene(shell_masks, tumor_mask, ct_node, patient_id, logfn=print):
    node_name = "Rims_%s" % patient_id
    old = slicer.mrmlScene.GetFirstNodeByName(node_name)
    if old:
        slicer.mrmlScene.RemoveNode(old)
    seg_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", node_name)
    seg_node.CreateDefaultDisplayNodes()
    seg_node.SetReferenceImageGeometryParameterFromVolumeNode(ct_node)

    def push(sitk_mask, label, colour):
        tmp = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLLabelMapVolumeNode", "__tmp_vis__")
        sitkUtils.PushVolumeToSlicer(sitk_mask, tmp)
        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(tmp, seg_node)
        seg  = seg_node.GetSegmentation()
        sid  = seg.GetNthSegmentID(seg.GetNumberOfSegments() - 1)
        segg = seg.GetSegment(sid)
        segg.SetName(label)
        segg.SetColor(*colour)
        slicer.mrmlScene.RemoveNode(tmp)

    push(tumor_mask, "Tumor (reference)", RIM_COLOURS["Tumor (reference)"])
    for name, _, _, colour in SHELLS:
        if name in shell_masks:
            push(shell_masks[name], name, RIM_COLOURS[name])

    seg_node.CreateClosedSurfaceRepresentation()
    slicer.util.setSliceViewerLayers(background=ct_node)
    slicer.util.resetSliceViews()
    slicer.util.resetThreeDViews()
    logfn("Rims shown as '%s'" % node_name)
    return seg_node


def show_habitats_in_scene(habitat_masks, ct_node, patient_id, logfn=print):
    node_name = "Habitats_%s" % patient_id
    old = slicer.mrmlScene.GetFirstNodeByName(node_name)
    if old:
        slicer.mrmlScene.RemoveNode(old)
    seg_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", node_name)
    seg_node.CreateDefaultDisplayNodes()
    seg_node.SetReferenceImageGeometryParameterFromVolumeNode(ct_node)
    colour_list = list(HABITAT_COLOURS.values())

    def push(sitk_mask, label, colour):
        tmp = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLLabelMapVolumeNode", "__tmp_hab__")
        sitkUtils.PushVolumeToSlicer(sitk_mask, tmp)
        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(tmp, seg_node)
        seg  = seg_node.GetSegmentation()
        sid  = seg.GetNthSegmentID(seg.GetNumberOfSegments() - 1)
        segg = seg.GetSegment(sid)
        segg.SetName(label)
        segg.SetColor(*colour)
        slicer.mrmlScene.RemoveNode(tmp)

    for i, (hname, hmask) in enumerate(sorted(habitat_masks.items())):
        push(hmask, hname, colour_list[i % len(colour_list)])

    seg_node.CreateClosedSurfaceRepresentation()
    slicer.util.resetThreeDViews()
    logfn("Habitats shown as '%s'" % node_name)
    return seg_node


def save_rims_to_file(shell_masks, tumor_mask, out_dir, patient_id, phase, logfn=print):
    os.makedirs(out_dir, exist_ok=True)
    safe_phase = phase.replace(" ", "_")
    ts    = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    fpath = os.path.join(out_dir, "%s_%s_%s_rims.seg.nrrd" % (patient_id, safe_phase, ts))
    tmp_seg = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "__tmp_save__")
    tmp_seg.CreateDefaultDisplayNodes()

    def push(sitk_mask, label, colour):
        tmp_lm = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLLabelMapVolumeNode", "__tmp_lm__")
        sitkUtils.PushVolumeToSlicer(sitk_mask, tmp_lm)
        slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(tmp_lm, tmp_seg)
        seg  = tmp_seg.GetSegmentation()
        sid  = seg.GetNthSegmentID(seg.GetNumberOfSegments() - 1)
        segg = seg.GetSegment(sid)
        segg.SetName(label)
        segg.SetColor(*colour)
        slicer.mrmlScene.RemoveNode(tmp_lm)

    push(tumor_mask, "Tumor (reference)", RIM_COLOURS["Tumor (reference)"])
    for name, _, _, colour in SHELLS:
        if name in shell_masks:
            push(shell_masks[name], name, RIM_COLOURS[name])
    slicer.util.saveNode(tmp_seg, fpath)
    slicer.mrmlScene.RemoveNode(tmp_seg)
    logfn("Rims saved: " + fpath)
    return fpath


# ======================================================================
#  SECTION 5 — PREPARE MASKS  (main thread only — ITK is not thread-safe)
# ======================================================================

def prepare_masks(ct_name, tumor_name, liver_name, liver_label,
                  shells_on, logfn=print):
    """All SimpleITK/ITK operations run here, on the main thread. Returns
    pure-NumPy arrays and pre-built sitk images safe to pass to a worker thread."""
    ct_node = slicer.mrmlScene.GetFirstNodeByName(ct_name)
    if ct_node is None:
        raise ValueError("CT node not found: " + ct_name)

    logfn("Pulling volumes from scene...")
    ct_itk    = sitk.Cast(sitkUtils.PullVolumeFromSlicer(ct_node), sitk.sitkFloat32)
    tumor_itk = pull_mask_as_sitk(tumor_name, ct_node)
    liver_raw = pull_mask_as_sitk(liver_name, ct_node)

    liver_itk = liver_raw if " > " in liver_name else \
                sitk.BinaryThreshold(liver_raw, liver_label, liver_label, 1, 0)

    logfn("Resampling to 1 mm isotropic (1-2 min)...")
    slicer.app.processEvents()
    ct    = resample_to_1mm(ct_itk)
    tmask = resample_mask(tumor_itk, ct)
    lmask = resample_mask(liver_itk, ct)

    t_arr = sitk.GetArrayFromImage(tmask)
    l_arr = sitk.GetArrayFromImage(lmask).astype(bool)
    tvox  = int(np.sum(t_arr))
    tvol  = tvox * float(np.prod(ct.GetSpacing()))
    lvox  = int(np.sum(l_arr))

    logfn("  Tumor : %d vox | %.0f mm3 | %.2f cc" % (tvox, tvol, tvol / 1000))
    logfn("  Liver : %d vox" % lvox)

    if tvox == 0:
        raise ValueError("Tumor mask empty after resampling.")
    if lvox < 100:
        raise ValueError("Liver mask empty. Try Liver label # = 1 or 2.")

    logfn("Normalising intensity (liver-mean subtraction)...")
    clamped = sitk.Clamp(ct, lowerBound=-100.0, upperBound=400.0)
    ct_arr  = sitk.GetArrayFromImage(clamped)
    lmean   = float(np.mean(ct_arr[l_arr]))
    logfn("  Liver mean HU (post-clamp): %.1f" % lmean)

    ct_norm = sitk.GetImageFromArray(ct_arr - lmean)
    ct_norm.CopyInformation(ct)

    logfn("Building peritumoral shells (single distance transform)...")
    slicer.app.processEvents()
    shell_masks = make_shells_fast(tmask, lmask, shells_on, logfn)

    if not shell_masks:
        raise ValueError("All shells empty. Tumor may be at the liver edge.")

    return ct_node, ct_norm, tmask, shell_masks, tvol, lmean


# ======================================================================
#  SECTION 6 — EXTRACTION  (heavy NumPy/radiomics work → worker thread)
# ======================================================================

def run_extraction_worker(ct_norm, tmask, shell_masks, tvol, lmean,
                           patient_id, phase, bin_width,
                           run_habitat, n_habitats, use_texture,
                           hab_fast, out_dir, save_rims_flag,
                           shell_masks_orig,
                           logfn=print):
    """
    Pure NumPy + PyRadiomics work — safe to run on a worker thread.
    No slicer/sitkUtils/Qt calls here.
    Returns (csv_path, feats, derived, habitat_masks, habitat_metrics).
    """
    ex_s = make_extractor(bin_width, shape=True)
    ex_n = make_extractor(bin_width, shape=False)
    feats = {}

    # ---- A. Intratumoral ----
    logfn("\nExtracting intratumoral features...")
    raw = ex_s.execute(ct_norm, tmask)
    feats["intra"] = {k: float(v) for k, v in raw.items() if k.startswith("original")}
    logfn("  %d features" % len(feats["intra"]))

    # ---- B. Peritumoral shells ----
    for name, mask in shell_masks.items():
        logfn("Extracting %s..." % name)
        raw = ex_n.execute(ct_norm, mask)
        feats[name] = {k: float(v) for k, v in raw.items() if k.startswith("original")}
        logfn("  %d features" % len(feats[name]))

    # ---- C. Rim derived features ----
    logfn("\nComputing derived rim features...")
    derived  = {}
    pfx_map  = {"peri_0-3mm": "0-3", "peri_3-5mm": "3-5", "peri_5-10mm": "5-10"}

    for zone, pfx in pfx_map.items():
        if zone not in feats:
            continue
        shared = set(feats["intra"]) & set(feats[zone])
        for k in shared:
            a, b = feats["intra"][k], feats[zone][k]
            s    = k.replace("original_", "")
            derived["diff_intra_vs_%s_%s"  % (pfx, s)] = a - b
            derived["ratio_intra_vs_%s_%s" % (pfx, s)] = (
                a / b if abs(b) > 1e-6 else float("nan"))

    zones = [z for z in ["peri_0-3mm", "peri_3-5mm", "peri_5-10mm"] if z in feats]
    for i in range(len(zones)):
        for j in range(i + 1, len(zones)):
            za, zb = zones[i], zones[j]
            pa, pb = pfx_map[za], pfx_map[zb]
            shared = set(feats[za]) & set(feats[zb])
            for k in shared:
                a, b = feats[za][k], feats[zb][k]
                s    = k.replace("original_", "")
                derived["diff_%s_vs_%s_%s"  % (pa, pb, s)] = a - b
                derived["ratio_%s_vs_%s_%s" % (pa, pb, s)] = (
                    a / b if abs(b) > 1e-6 else float("nan"))

    m03  = feats.get("peri_0-3mm",  {}).get("original_firstorder_Mean")
    m510 = feats.get("peri_5-10mm", {}).get("original_firstorder_Mean")
    if m03 is not None and m510 is not None:
        derived["gradient_infiltration_index"] = m03 - m510
        logfn("  GII : %.4f" % derived["gradient_infiltration_index"])

    im_norm = feats.get("intra", {}).get("original_firstorder_Mean")
    pm_norm = feats.get("peri_0-3mm", {}).get("original_firstorder_Mean")
    if im_norm is not None and pm_norm is not None:
        im_abs = im_norm + lmean
        if abs(im_abs) > 1.0:
            derived["peripheral_enhancement_ratio_0-3mm"] = (pm_norm + lmean) / im_abs
            logfn("  PER : %.4f" % derived["peripheral_enhancement_ratio_0-3mm"])

    ie = feats.get("intra", {}).get("original_firstorder_Energy")
    re_ = feats.get("peri_0-3mm", {}).get("original_firstorder_Energy")
    if ie and re_ and (ie + re_) > 0:
        derived["rim_energy_fraction_0-3mm"] = re_ / (ie + re_)

    # ---- D. Habitat Imaging ----
    habitat_masks   = {}
    habitat_metrics = {}

    if run_habitat:
        logfn("\n--- Habitat Imaging (k=%d, texture=%s) ---" % (n_habitats, use_texture))
        try:
            habitat_masks, cluster_stats = compute_habitats(
                ct_norm, tmask, n_habitats=n_habitats,
                use_texture=use_texture, logfn=logfn)

            ex_hab   = make_extractor(bin_width, shape=False, fast=hab_fast)
            mode_str = "1st-order only" if hab_fast else "full texture"
            logfn("  Feature mode: %s" % mode_str)
            for hname, hmask in sorted(habitat_masks.items()):
                hvox = int(np.sum(sitk.GetArrayFromImage(hmask)))
                if hvox < 50:
                    logfn("  %s: SKIPPED (too small)" % hname)
                    continue
                logfn("  Extracting %s..." % hname)
                raw = ex_hab.execute(ct_norm, hmask)
                feats[hname] = {k: float(v) for k, v in raw.items()
                                if k.startswith("original")}
                logfn("    %d features" % len(feats[hname]))

            habitat_metrics = compute_habitat_metrics(
                cluster_stats, n_habitats, ct_norm, tmask)
            logfn("\n  Habitat metrics:")
            for k, v in habitat_metrics.items():
                logfn("    %s : %.4f" % (k, v))

            # Cross-habitat derived features
            h_keys = sorted([k for k in feats if k.startswith("habitat_")])
            if len(h_keys) >= 2:
                h_low, h_high = h_keys[0], h_keys[-1]
                shared = set(feats.get(h_low, {}).keys()) & set(feats.get(h_high, {}).keys())
                for k in shared:
                    s = k.replace("original_", "")
                    a, b = feats[h_low][k], feats[h_high][k]
                    derived["habitat_lowhigh_diff_%s"  % s] = a - b
                    derived["habitat_lowhigh_ratio_%s" % s] = (
                        a / b if abs(b) > 1e-6 else float("nan"))

        except Exception as e:
            logfn("WARNING: Habitat imaging failed — %s" % str(e))
            logfn(traceback.format_exc())
            logfn("Continuing with rim features only.")

    # ---- E. Save CSV ----
    os.makedirs(out_dir, exist_ok=True)
    ts         = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_phase = phase.replace(" ", "_")
    csv_path   = os.path.join(out_dir, "%s_%s_%s.csv" % (patient_id, safe_phase, ts))

    row = {"patient_id": patient_id, "phase": phase,
           "tumor_vol_mm3": tvol, "liver_mean_hu": lmean}
    for zone, fd in feats.items():
        for k, v in fd.items():
            row["%s_%s" % (zone, k)] = v
    row.update(derived)
    row.update({"hab_" + k: v for k, v in habitat_metrics.items()})
    pd.DataFrame([row]).to_csv(csv_path, index=False)

    logfn("\n" + "=" * 54)
    logfn("CSV DONE — %d total features" % len(row))
    logfn(csv_path)
    logfn("=" * 54)

    return csv_path, shell_masks, habitat_masks


# ======================================================================
#  SECTION 7 — GUI
# ======================================================================

def launch():
    vol_names  = get_volume_names()
    mask_names = get_mask_names()

    if not vol_names and not mask_names:
        slicer.util.warningDisplay(
            "No volumes found in the scene.\n\n"
            "Load your CT, tumor mask, and liver mask first via\n"
            "File -> Add Data, then re-run the script.",
            windowTitle="Radiomics Extraction")
        return

    dlg = qt.QDialog(slicer.util.mainWindow())
    dlg.setWindowTitle("Peritumoral Radiomics + Habitat Imaging")
    dlg.setMinimumWidth(590)
    dlg.resize(620, 840)
    dlg.setModal(False)
    dlg.setWindowFlags(
        qt.Qt.Window |
        qt.Qt.WindowTitleHint |
        qt.Qt.WindowSystemMenuHint |
        qt.Qt.WindowMinimizeButtonHint |
        qt.Qt.WindowMaximizeButtonHint |
        qt.Qt.WindowCloseButtonHint)

    outer = qt.QVBoxLayout(dlg)
    outer.setSpacing(0)
    outer.setContentsMargins(0, 0, 0, 0)

    scroll = qt.QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(qt.QFrame.NoFrame)
    scroll.setHorizontalScrollBarPolicy(qt.Qt.ScrollBarAlwaysOff)
    inner_widget = qt.QWidget()
    main = qt.QVBoxLayout(inner_widget)
    main.setSpacing(8)
    main.setContentsMargins(16, 14, 16, 8)
    scroll.setWidget(inner_widget)
    outer.addWidget(scroll, 1)

    def section(txt):
        lb = qt.QLabel("  " + txt)
        lb.setStyleSheet(
            "background:#1a5276;color:white;font-weight:bold;"
            "font-size:11px;padding:3px 6px;border-radius:3px;")
        main.addWidget(lb)

    def add_combo_row(label_txt, items, tip=""):
        row = qt.QHBoxLayout()
        lbl = qt.QLabel(label_txt)
        lbl.setFixedWidth(155)
        lbl.setStyleSheet("font-weight:bold;font-size:12px;")
        cb = qt.QComboBox()
        cb.setMinimumWidth(340)
        for item in items:
            cb.addItem(item)
        if tip:
            cb.setToolTip(tip)
        row.addWidget(lbl)
        row.addWidget(cb)
        main.addLayout(row)
        return cb

    def styled_btn(txt, color, w=None, h=None):
        b = qt.QPushButton(txt)
        if w: b.setFixedWidth(w)
        if h: b.setFixedHeight(h)
        b.setStyleSheet(
            "QPushButton{background:" + color + ";color:white;border-radius:3px;"
            "font-weight:bold;padding:5px 10px;}"
            "QPushButton:hover{background:#1a5276;}"
            "QPushButton:disabled{background:#aaa;}")
        return b

    # Title
    ttl = qt.QLabel("Peritumoral Radiomics + Habitat Imaging")
    ttl.setStyleSheet("font-size:15px;font-weight:bold;color:#1a5276;")
    main.addWidget(ttl)
    sub = qt.QLabel("Rim features  +  intra-tumoral habitat partitioning  +  combined CSV")
    sub.setStyleSheet("color:#666;font-size:10px;margin-bottom:4px;")
    main.addWidget(sub)

    # ---- Volume selection ----
    section("Select Volumes from Scene")
    ct_cb    = add_combo_row("CT Volume:",  vol_names)
    tumor_cb = add_combo_row("Tumor Mask:", mask_names)
    liver_cb = add_combo_row("Liver Mask:", mask_names)

    for i, n in enumerate(vol_names):
        if any(w in n.lower() for w in ["scan", "ct", "art", "portal", "phase", "series"]):
            ct_cb.setCurrentIndex(i); break
    for i, n in enumerate(mask_names):
        if any(w in n.lower() for w in ["tumor", "tumour", "lesion", "mass"]):
            tumor_cb.setCurrentIndex(i); break
    for i, n in enumerate(mask_names):
        if any(w in n.lower() for w in ["liver", "hepat"]):
            liver_cb.setCurrentIndex(i); break

    lrow = qt.QHBoxLayout()
    llbl = qt.QLabel("Liver label #:")
    llbl.setFixedWidth(155)
    llbl.setStyleSheet("font-weight:bold;font-size:12px;")
    liver_spin = qt.QSpinBox()
    liver_spin.setRange(1, 20)
    liver_spin.setValue(1)
    liver_spin.setFixedWidth(58)
    liver_tip = qt.QLabel("")
    liver_tip.setStyleSheet("font-size:10px;font-style:italic;")
    lrow.addWidget(llbl)
    lrow.addWidget(liver_spin)
    lrow.addWidget(liver_tip)
    lrow.addStretch()
    main.addLayout(lrow)

    def update_liver_tip():
        if " > " in liver_cb.currentText:
            liver_spin.setEnabled(False)
            liver_tip.setText("Not needed — segment already binary  ✓")
            liver_tip.setStyleSheet("color:#27ae60;font-size:10px;")
        else:
            liver_spin.setEnabled(True)
            liver_tip.setText("Label index in LabelMap (usually 1)")
            liver_tip.setStyleSheet("color:#888;font-size:10px;font-style:italic;")
    liver_cb.connect("currentIndexChanged(int)", lambda _=None: update_liver_tip())
    update_liver_tip()

    # ---- Rim zones ----
    section("Rim Zones")
    rim_row = qt.QHBoxLayout()
    cb03    = qt.QCheckBox("0–3 mm");  cb03.setChecked(True)
    cb35    = qt.QCheckBox("3–5 mm");  cb35.setChecked(True)
    cb510   = qt.QCheckBox("5–10 mm"); cb510.setChecked(True)
    btn_preview = styled_btn("Preview Rims", "#6c3483", w=135)
    for w in [cb03, cb35, cb510, btn_preview]:
        rim_row.addWidget(w)
    rim_row.addStretch()
    main.addLayout(rim_row)
    legend = qt.QLabel("  Blue=Tumor  Green=0–3mm  Orange=3–5mm  Red=5–10mm")
    legend.setStyleSheet("color:#555;font-size:10px;font-style:italic;")
    main.addWidget(legend)

    # ---- Habitat Imaging ----
    section("Habitat Imaging  (intra-tumoral subregion partitioning)")

    hab_en_row = qt.QHBoxLayout()
    cb_habitat = qt.QCheckBox("Enable habitat imaging")
    cb_habitat.setChecked(True)
    cb_habitat.setStyleSheet("font-weight:bold;font-size:12px;")
    hab_en_row.addWidget(cb_habitat)
    hab_en_row.addStretch()
    main.addLayout(hab_en_row)

    hab_settings = qt.QWidget()
    hab_layout   = qt.QVBoxLayout(hab_settings)
    hab_layout.setContentsMargins(0, 0, 0, 0)
    hab_layout.setSpacing(4)

    nrow = qt.QHBoxLayout()
    nlbl = qt.QLabel("Number of habitats (k):")
    nlbl.setFixedWidth(200)
    nlbl.setStyleSheet("font-size:12px;")
    n_hab_spin = qt.QSpinBox()
    n_hab_spin.setRange(2, 4)
    n_hab_spin.setValue(3)
    n_hab_spin.setFixedWidth(58)
    n_hab_spin.setToolTip("k=3 recommended: low=necrotic, mid=viable, high=enhancing")
    nrow.addWidget(nlbl)
    nrow.addWidget(n_hab_spin)
    nrow.addStretch()
    hab_layout.addLayout(nrow)

    tex_row    = qt.QHBoxLayout()
    cb_texture = qt.QCheckBox("Enrich with local texture (3×3×3 std + gradient)  [fast]")
    cb_texture.setChecked(True)
    cb_texture.setStyleSheet("font-size:11px;")
    cb_texture.setToolTip(
        "Uses fast uniform_filter std (not generic_filter) — adds ~3 seconds.\n"
        "Makes habitats correspond to tissue biology rather than bulk HU alone.")
    tex_row.addWidget(cb_texture)
    tex_row.addStretch()
    hab_layout.addLayout(tex_row)

    show_hab_row = qt.QHBoxLayout()
    cb_show_hab  = qt.QCheckBox("Show habitat segments in 3D viewer after extraction")
    cb_show_hab.setChecked(True)
    cb_show_hab.setStyleSheet("font-size:11px;")
    show_hab_row.addWidget(cb_show_hab)
    show_hab_row.addStretch()
    hab_layout.addLayout(show_hab_row)

    fast_row    = qt.QHBoxLayout()
    cb_hab_fast = qt.QCheckBox("Fast mode: first-order features only for habitats  (recommended)")
    cb_hab_fast.setChecked(True)
    cb_hab_fast.setStyleSheet("font-size:11px;")
    cb_hab_fast.setToolTip(
        "ON (default): first-order only — ~30 sec total.\n"
        "OFF: full GLCM/GLSZM per habitat — 10-20 min, can still be slow.")
    fast_row.addWidget(cb_hab_fast)
    fast_row.addStretch()
    hab_layout.addLayout(fast_row)

    hab_hint = qt.QLabel(
        "  k=3 colours:  Blue=low HU (necrosis)  Yellow=mid (viable)  Red=high (enhancing)")
    hab_hint.setStyleSheet("color:#555;font-size:10px;font-style:italic;")
    hab_layout.addWidget(hab_hint)
    main.addWidget(hab_settings)

    def toggle_habitat_settings():
        hab_settings.setEnabled(cb_habitat.isChecked())
    cb_habitat.connect("toggled(bool)", lambda _=None: toggle_habitat_settings())
    toggle_habitat_settings()

    # ---- Patient info ----
    section("Patient Info")
    prow = qt.QHBoxLayout()
    plbl = qt.QLabel("Patient ID:")
    plbl.setFixedWidth(155)
    plbl.setStyleSheet("font-weight:bold;font-size:12px;")
    pid_ed = qt.QLineEdit()
    pid_ed.setPlaceholderText("e.g. 35")
    pid_ed.setFixedWidth(200)
    prow.addWidget(plbl)
    prow.addWidget(pid_ed)
    prow.addStretch()
    main.addLayout(prow)

    crow = qt.QHBoxLayout()
    clbl = qt.QLabel("CT Phase:")
    clbl.setFixedWidth(155)
    clbl.setStyleSheet("font-weight:bold;font-size:12px;")
    phase_cb = qt.QComboBox()
    for p in ["Arterial", "Portal_Venous", "Delayed", "Non_contrast", "Equilibrium"]:
        phase_cb.addItem(p)
    phase_cb.setFixedWidth(200)
    crow.addWidget(clbl)
    crow.addWidget(phase_cb)
    crow.addStretch()
    main.addLayout(crow)

    # ---- Settings ----
    section("Settings")
    srow = qt.QHBoxLayout()
    slbl = qt.QLabel("Bin width (HU):")
    slbl.setFixedWidth(155)
    slbl.setStyleSheet("font-weight:bold;font-size:12px;")
    bin_sb = qt.QSpinBox()
    bin_sb.setRange(5, 100)
    bin_sb.setValue(25)
    bin_sb.setFixedWidth(58)
    srow.addWidget(slbl)
    srow.addWidget(bin_sb)
    srow.addStretch()
    main.addLayout(srow)

    # ---- Output ----
    section("Output")
    orow = qt.QHBoxLayout()
    olbl = qt.QLabel("Save folder:")
    olbl.setFixedWidth(155)
    olbl.setStyleSheet("font-weight:bold;font-size:12px;")
    out_ed  = qt.QLineEdit()
    out_ed.setText("")
    out_btn = styled_btn("Browse...", "#2471a3", w=75)

    def on_browse(_=None):
        path = qt.QFileDialog.getExistingDirectory(dlg, "Select output folder")
        if isinstance(path, (list, tuple)):
            path = path[0]
        if path:
            out_ed.setText(path)

    out_btn.connect("clicked(bool)", on_browse)
    orow.addWidget(olbl)
    orow.addWidget(out_ed)
    orow.addWidget(out_btn)
    main.addLayout(orow)

    sav_row      = qt.QHBoxLayout()
    cb_show      = qt.QCheckBox("Show rims in viewer after extraction")
    cb_save_rims = qt.QCheckBox("Save rims as .seg.nrrd")
    cb_show.setChecked(True)
    cb_save_rims.setChecked(True)
    for w in [cb_show, cb_save_rims]:
        w.setStyleSheet("font-size:11px;")
    sav_row.addWidget(cb_show)
    sav_row.addSpacing(18)
    sav_row.addWidget(cb_save_rims)
    sav_row.addStretch()
    main.addLayout(sav_row)

    # ---- Log ----
    log_header = qt.QLabel("  Log")
    log_header.setStyleSheet(
        "background:#1a5276;color:white;font-weight:bold;"
        "font-size:11px;padding:3px 6px;")
    outer.addWidget(log_header)

    log_box = qt.QPlainTextEdit()
    log_box.setReadOnly(True)
    log_box.setFixedHeight(160)
    log_box.setStyleSheet(
        "background:#111827;color:#d1e8ff;"
        "font-family:Consolas,Courier,monospace;"
        "font-size:11px;border:none;padding:5px;")
    outer.addWidget(log_box)

    def W(msg):
        log_box.appendPlainText(str(msg))
        slicer.app.processEvents()

    # ---- Bottom buttons ----
    bot_widget = qt.QWidget()
    bot_widget.setStyleSheet("background:#f0f0f0;border-top:1px solid #ccc;")
    bot = qt.QHBoxLayout(bot_widget)
    bot.setContentsMargins(12, 8, 12, 10)
    btn_refresh  = styled_btn("Refresh",       "#5d6d7e", w=80)
    btn_save_now = styled_btn("Save Rims Now", "#8e44ad", w=120)
    btn_save_now.setEnabled(False)
    btn_open     = styled_btn("Open Folder",   "#27ae60", w=105)
    btn_open.setEnabled(False)
    btn_run = qt.QPushButton("▶  Run Extraction")
    btn_run.setFixedHeight(38)
    btn_run.setStyleSheet(
        "QPushButton{background:#c0392b;color:white;font-size:14px;"
        "font-weight:bold;border-radius:4px;padding:5px 20px;}"
        "QPushButton:hover{background:#a93226;}"
        "QPushButton:disabled{background:#999;}")
    for w in [btn_refresh, btn_save_now, btn_open]:
        bot.addWidget(w)
    bot.addStretch()
    bot.addWidget(btn_run)
    outer.addWidget(bot_widget)

    _state = {
        "last_out":      None,
        "shell_masks":   None,
        "habitat_masks": None,
        "tumor_mask":    None,
        "ct_node":       None,
    }

    def get_shells_on():
        return (
            (["peri_0-3mm"]  if cb03.isChecked()  else []) +
            (["peri_3-5mm"]  if cb35.isChecked()  else []) +
            (["peri_5-10mm"] if cb510.isChecked() else [])
        )

    def validate(need_output=True):
        errs = []
        if not ct_cb.currentText:    errs.append("Select a CT volume.")
        if not tumor_cb.currentText: errs.append("Select a tumor mask.")
        if not liver_cb.currentText: errs.append("Select a liver mask.")
        if need_output and not out_ed.text.strip():
            errs.append("Set an output folder.")
        if not (cb03.isChecked() or cb35.isChecked() or cb510.isChecked()):
            errs.append("Tick at least one rim zone.")
        return errs

    def on_refresh(_=None):
        nv = get_volume_names()
        nm = get_mask_names()
        for prev, cb, items in [
            (ct_cb.currentText,    ct_cb,    nv),
            (tumor_cb.currentText, tumor_cb, nm),
            (liver_cb.currentText, liver_cb, nm),
        ]:
            cb.clear()
            for x in items: cb.addItem(x)
            idx = cb.findText(prev)
            if idx >= 0: cb.setCurrentIndex(idx)
        update_liver_tip()
        W("Refreshed — %d CT volumes, %d mask entries" % (len(nv), len(nm)))

    btn_refresh.connect("clicked(bool)", on_refresh)

    def on_preview(_=None):
        errs = validate(need_output=False)
        if errs:
            W("ERROR: " + "  |  ".join(errs)); return
        btn_preview.setEnabled(False)
        btn_preview.setText("Building...")
        slicer.app.processEvents()
        W("\n--- Preview Rims ---")
        try:
            ct_node, ct_norm, tmask, shell_masks, tvol, lmean = \
                prepare_masks(ct_cb.currentText, tumor_cb.currentText,
                              liver_cb.currentText, liver_spin.value,
                              get_shells_on(), W)
            _state["shell_masks"] = shell_masks
            _state["tumor_mask"]  = tmask
            _state["ct_node"]     = ct_node
            show_rims_in_scene(shell_masks, tmask, ct_node,
                               pid_ed.text.strip() or "patient", W)
            btn_save_now.setEnabled(True)
            btn_preview.setText("Rims Shown ✓")
            btn_preview.setStyleSheet(
                "QPushButton{background:#27ae60;color:white;border-radius:3px;"
                "font-weight:bold;padding:5px 10px;}"
                "QPushButton:hover{background:#1e8449;}")
        except Exception as e:
            W("ERROR: " + str(e))
            W(traceback.format_exc())
            btn_preview.setText("Preview Rims")
            btn_preview.setStyleSheet(
                "QPushButton{background:#6c3483;color:white;border-radius:3px;"
                "font-weight:bold;padding:5px 10px;}"
                "QPushButton:hover{background:#1a5276;}")
        finally:
            btn_preview.setEnabled(True)

    btn_preview.connect("clicked(bool)", on_preview)

    def on_save_rims_now(_=None):
        if not _state["shell_masks"]:
            W("ERROR: No rims built yet. Click Preview Rims first."); return
        out = out_ed.text.strip()
        if not out:
            W("ERROR: Set an output folder first."); return
        try:
            p = save_rims_to_file(_state["shell_masks"], _state["tumor_mask"],
                                  out, pid_ed.text.strip() or "patient",
                                  phase_cb.currentText, W)
            _state["last_out"] = os.path.dirname(p)
            btn_open.setEnabled(True)
        except Exception as e:
            W("ERROR: " + str(e))

    btn_save_now.connect("clicked(bool)", on_save_rims_now)

    def on_open(_=None):
        if _state["last_out"]:
            qt.QDesktopServices.openUrl(qt.QUrl.fromLocalFile(_state["last_out"]))
    btn_open.connect("clicked(bool)", on_open)

    def on_run(_=None):
        errs = validate()
        if errs:
            W("ERROR: " + "  |  ".join(errs)); return

        btn_run.setEnabled(False)
        btn_run.setText("Running…")
        btn_open.setEnabled(False)
        W("\n" + "=" * 54)
        W("Step 1/2: Preparing masks (main thread — ITK)...")
        slicer.app.processEvents()

        # ── STEP 1: All ITK work on main thread ────────────────────────
        try:
            ct_node, ct_norm, tmask, shell_masks, tvol, lmean = \
                prepare_masks(ct_cb.currentText, tumor_cb.currentText,
                              liver_cb.currentText, liver_spin.value,
                              [s for s in ["peri_0-3mm", "peri_3-5mm", "peri_5-10mm"]
                               if (s == "peri_0-3mm" and cb03.isChecked()) or
                                  (s == "peri_3-5mm" and cb35.isChecked()) or
                                  (s == "peri_5-10mm" and cb510.isChecked())],
                              W)
        except Exception as e:
            W("ERROR in mask preparation: " + str(e))
            W(traceback.format_exc())
            btn_run.setEnabled(True)
            btn_run.setText("▶  Run Extraction")
            return

        W("Step 2/2: Extracting features (worker thread — NumPy/radiomics)...")
        slicer.app.processEvents()

        # Snapshot GUI values before handing off to worker
        params = dict(
            ct_norm      = ct_norm,
            tmask        = tmask,
            shell_masks  = shell_masks,
            tvol         = tvol,
            lmean        = lmean,
            patient_id   = pid_ed.text.strip() or "patient",
            phase        = phase_cb.currentText,
            bin_width    = bin_sb.value,
            run_habitat  = cb_habitat.isChecked(),
            n_habitats   = n_hab_spin.value,
            use_texture  = cb_texture.isChecked(),
            hab_fast     = cb_hab_fast.isChecked(),
            out_dir      = out_ed.text.strip(),
            save_rims_flag = cb_save_rims.isChecked(),
            shell_masks_orig = shell_masks,
        )

        # ── STEP 2: Pure NumPy/sklearn/radiomics → worker thread ──────────
        log_q   = _queue.Queue()
        result_b = [None]
        error_b  = [None]
        done_ev  = threading.Event()

        def threadsafe_log(msg):
            log_q.put(str(msg))

        def worker():
            try:
                result_b[0] = run_extraction_worker(logfn=threadsafe_log, **params)
            except Exception as exc:
                error_b[0] = exc
            finally:
                done_ev.set()

        threading.Thread(target=worker, daemon=True).start()

        # Pump UI while worker runs
        while not done_ev.wait(timeout=0.08):
            while not log_q.empty():
                W(log_q.get_nowait())
            slicer.app.processEvents()
        while not log_q.empty():
            W(log_q.get_nowait())
        slicer.app.processEvents()

        # ── STEP 3: Visualisation back on main thread ─────────────────────
        if error_b[0] is not None:
            W("ERROR: " + str(error_b[0]))
            W(traceback.format_exc())
            btn_run.setText("Error — see log")
            btn_run.setStyleSheet(
                "QPushButton{background:#c0392b;color:white;font-size:14px;"
                "font-weight:bold;border-radius:4px;padding:5px 20px;}")
        else:
            csv_path, shell_masks_out, habitat_masks = result_b[0]
            _state.update({
                "last_out":      os.path.dirname(csv_path),
                "shell_masks":   shell_masks_out,
                "habitat_masks": habitat_masks,
                "tumor_mask":    tmask,
                "ct_node":       ct_node,
            })
            if cb_show.isChecked() and shell_masks_out:
                W("Loading rims into viewer...")
                show_rims_in_scene(shell_masks_out, tmask, ct_node,
                                   params["patient_id"], W)
            if cb_show_hab.isChecked() and cb_habitat.isChecked() and habitat_masks:
                W("Loading habitats into viewer...")
                show_habitats_in_scene(habitat_masks, ct_node, params["patient_id"], W)
            if cb_save_rims.isChecked() and shell_masks_out:
                save_rims_to_file(shell_masks_out, tmask, params["out_dir"],
                                  params["patient_id"], params["phase"], W)
            btn_save_now.setEnabled(True)
            btn_open.setEnabled(True)
            btn_run.setText("Done ✓")
            btn_run.setStyleSheet(
                "QPushButton{background:#27ae60;color:white;font-size:14px;"
                "font-weight:bold;border-radius:4px;padding:5px 20px;}"
                "QPushButton:hover{background:#1e8449;}")

        btn_run.setEnabled(True)
        slicer.app.processEvents()

    btn_run.connect("clicked(bool)", on_run)

    W("Scene — CT: %d   Masks: %d" % (len(vol_names), len(mask_names)))
    W("-" * 52)
    W("Workflow:")
    W("  1. Verify CT + Tumor + Liver dropdowns")
    W("  2. Optional: Preview Rims to inspect zones on CT")
    W("  3. Enable Habitat Imaging if needed (k=3 recommended)")
    W("  4. Click  Run Extraction")
    W("  NOTE: resampling runs on the main thread (~1-2 min),")
    W("        then radiomics/k-means run in background.")
    W("        Slicer stays responsive throughout.")

    dlg.show()
    return dlg


gui_window = launch()
