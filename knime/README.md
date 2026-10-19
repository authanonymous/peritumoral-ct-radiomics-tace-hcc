# KNIME workflow

`knime_workflow.knwf` (exported without data; KNIME Analytics Platform 5.8.3): reads the merged PyRadiomics feature table, applies preprocessing (type casting, joins against the clinical table, missing-value handling, duplicate/constant-column removal), the Low Variance Filter (threshold 0.09) and the Correlation Filter node (default settings, |r| > 0.90), and exports the 624-feature table.

Column counts at each step are in `../results/preprocessing_step_counts.csv`.

This copy has been sanitised for public release: the local file paths in the CSV Reader/Writer and Excel Reader nodes, and the Windows author/editor metadata in the workflow file, were replaced with generic placeholders. No node, setting, or data-processing logic was changed.
