# unsupflowhelper 1.1

## Changes
- Fixed a bug where more than one split.by argument for plot_density_umap gave an error
- Removed forced statistical comparisons on cluster_jitterplot
- Fixed scaling on filtered cluster heatmaps, so that colors remained identical to their unfiltered version
- Fixed several bugs for trajectory analysis functions
- Added FindTrajectory markers function
- Added pseudotime imputation function.
- Added trajectory analysis vignette
- Added RShiny GUI for selecting files and adding metadata. 

## Known bugs
- Issue on Windows OS where reading CSV files as input gives an error when trying to convert to a flowSet
- Issue when FCS files are empty.

# unsupflowhelpder 0.1

*This is the first release of unsupflowhelper
