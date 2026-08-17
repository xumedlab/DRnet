# Provenance of the Voigt 2019 author cluster mapping

## Source

- Article: Voigt AP et al. *Molecular characterization of foveal versus peripheral human retina by single-cell RNA sequencing*. Experimental Eye Research. 2019;184:234–242.
- DOI: `10.1016/j.exer.2019.05.001`
- Mapping location: Figure 1F on page 38 of the 42-page accepted-manuscript PDF.
- Verified: 13 August 2026. The source PDF is not redistributed because of copyright restrictions.

## Mapping evidence

The `author_label` field in `voigt2019_author_cluster_mapping.csv` was transcribed directly from Figure 1F:

- clusters 1–2: rods
- clusters 3–4: cones
- clusters 5–6: bipolar cells
- cluster 7: retinal ganglion cells
- cluster 8A: horizontal cells
- cluster 8B: amacrine cells
- cluster 9: unknown
- cluster 10: pericytes
- cluster 11: endothelial cells
- cluster 12: microglia
- clusters 13–17: glial cells

Results section 3.2 provides an independent textual cross-check. Cluster 9 lacked the selected cell-specific genes and was described as unknown; cluster 10 represented pericytes or smooth-muscle-like mural cells; cluster 11 represented endothelial cells; cluster 12 represented microglia; and clusters 13–17 were glial cells described as Müller cells and/or astrocytes.

## Analysis-label boundaries

- Cluster 9 is excluded from localization summaries and is not reclassified.
- The original label for clusters 13–17 is retained as `Glial cells` in the author-mapping table.
- The analysis label `Müller-enriched glia` reflects the source report that all five clusters expressed high `RLBP1` with relatively low `ALDH1L1` and `GFAP`. It is a source-constrained analysis label, not a claim that the authors assigned every cell specifically as a Müller cell.
- No marker-based reclustering, automated annotation, or inferred remapping of the authors' clusters was performed.
