@ECHO OFF
cmd /k "cd %~dp0 & env\scripts\activate & cd %~dp0 & python Mass_VisualizationChr1.py & rscript Plotting_Masses.R & deactivate
pause

