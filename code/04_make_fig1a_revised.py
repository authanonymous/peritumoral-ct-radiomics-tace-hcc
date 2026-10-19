import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import font_manager as fm
plt.rcParams['font.family']='DejaVu Serif'
W,H=1430,1324
fig=plt.figure(figsize=(W/100,H/100),dpi=100); ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,W); ax.set_ylim(H,0); ax.axis('off')
BL='#154F73'; BF='#EBF2FA'; PF='#FDECEA'; PB='#B03A2E'; PT='#922B21'; GR='#1E7B34'; GF='#E8F6EC'; AR='#2C3E50'
def box(x0,y0,x1,y1,fc,ec,lw=2.2): ax.add_patch(FancyBboxPatch((x0,y0),x1-x0,y1-y0,boxstyle='round,pad=0,rounding_size=14',fc=fc,ec=ec,lw=lw))
def txt(x,y,s,size=13,w='normal',st='normal',c='#222'): ax.text(x,y,s,fontsize=size,fontweight=w,fontstyle=st,color=c,ha='center',va='center')
LX0,LX1,RX0,RX1=35,930,966,1405; CX=(LX0+LX1)/2; RC=(RX0+RX1)/2
def arrow(p,q): ax.add_patch(FancyArrowPatch(p,q,arrowstyle='-|>',mutation_scale=22,lw=2.2,color=AR,shrinkA=0,shrinkB=0))
G='#5D6D7E'
box(LX0,35,LX1,165,BF,BL); txt(CX,70,'WAW-TACE Public Dataset',16,'bold',c=BL); txt(CX,104,'377 segmented HCC lesions \u00b7 233 patients',13); txt(CX,138,'Retrospective \u00b7 Pre-treatment CT \u00b7 May 2016\u2013April 2021 \u00b7 Publicly available (no IRB required)',10.5,st='italic',c=G)
box(LX0,305,LX1,450,BF,BL); txt(CX,334,'375 Lesions \u00b7 231 Patients',16,'bold',c=BL); txt(CX,364,'1,526 radiomic features extracted from all lesions',13); txt(CX,392,'PyRadiomics v3.0 \u00b7 3D Slicer v5.10.0 \u00b7 IBSI-compliant',10.5,st='italic',c=G); txt(CX,414,'Regions: intratumoral \u00b7 0\u20133 mm \u00b7 3\u20135 mm \u00b7 5\u201310 mm \u00b7 ratio/difference \u00b7 habitat (k-means, k=3)',10,st='italic',c=G); txt(CX,435,'Categories: first-order \u00b7 GLCM \u00b7 GLSZM \u00b7 GLRLM \u00b7 GLDM \u00b7 Shape (intratumoral only)',10,st='italic',c=G)
box(LX0,610,LX1,730,BF,BL); txt(CX,640,'231 Patients \u2014 Index Lesion Retained',16,'bold',c=BL); txt(CX,672,'Largest lesion per patient selected for patient-level survival analysis',13); txt(CX,703,'Peritumoral zone generation initiated \u00b7 3 non-overlapping shells: 0\u20133, 3\u20135, 5\u201310 mm',10,st='italic',c=G)
box(LX0,880,LX1,975,BF,BL); txt(CX,915,'213 Patients',16,'bold',c=BL); txt(CX,950,'Non-overlapping peritumoral zones generated (0\u20133, 3\u20135, 5\u201310 mm)',13)
box(LX0,1130,LX1,1290,GF,GR,2.6); ax.add_patch(FancyBboxPatch((LX0,1130),LX1-LX0,52,boxstyle='round,pad=0,rounding_size=14',fc=GR,ec=GR,lw=0)); txt(CX,1157,'Final Analytic Cohort \u2014 110 Patients',16,'bold',c='white')
txt(CX,1208,'85 men (77.3%) \u00b7 Mean age 65.9 \u00b1 9.9 years',13,'bold',c='#145A32'); txt(CX,1234,'75 events (68.2%) \u00b7 Median overall survival: 32.5 months',12.5,c='#1E8449'); txt(CX,1258,'BCLC A: 72.7% \u00b7 BCLC B: 27.3% \u00b7 Multifocal disease: 35.5% \u00b7 Child-Pugh category 1 (presumed class A): 85.5%',10,c='#1E8449'); txt(CX,1279,'Arterial-phase CT \u00b7 Index lesion per patient \u00b7 WAW-TACE dataset',10,st='italic',c=G)
def excl(j,y0,y1,l1,l2,l3=None):
    ax.plot([CX],[j],'o',color=AR,ms=7); arrow((CX,j),(RX0-4,j)); box(RX0,y0,RX1,y1,PF,PB)
    cy=(y0+y1)/2; n=2+(l3 is not None); ys=[cy-(n-1)*14+28*k for k in range(n)]
    txt(RC,ys[0],l1,13,'bold',c=PT); txt(RC,ys[1],l2,11,c=PT)
    if l3: txt(RC,ys[2],l3,11,c=PT)
arrow((CX,165),(CX,303)); excl(245,185,305,'Excluded: 2 lesions from 2 patients','Non-contrast CT segmentations only','(no arterial-phase data available)')
arrow((CX,450),(CX,608)); excl(540,480,600,'Excluded: 144 non-index lesions','Multifocal disease: single largest (index) lesion retained','per patient for patient-level survival modelling')
arrow((CX,730),(CX,878)); excl(805,740,870,'Excluded: 18 patients','Peritumoral rims overlapping the shell or tumour mask','of a neighbouring lesion')
arrow((CX,975),(CX,1128)); excl(1050,990,1110,'Excluded: 103 patients','Non-arterial-phase CT segmentations','(portal venous or other imaging phases)')
txt(RC,85,'EXCLUSION CRITERIA',12,'bold',c=PT); ax.plot([RX0,RX1],[103,103],color=PB,lw=1.5,ls=(0,(4,3)))
fig.savefig('../results/Fig1A_revised.png',dpi=100)
from PIL import Image
im=Image.open('../results/Fig1A_revised.png').convert('RGB'); print(im.size)
im.save('../results/Fig1A_revised.tif',compression='tiff_lzw'); im.convert('L').save('../results/Fig1A_revised_BW.png')
