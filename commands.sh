sbatch toJobMetricProcess.sh valencia day mean-max historical pr True
sbatch toJobMetricProcess.sh valencia day mean-max ssp370 pr True
sbatch toJobMetricProcess.sh valencia day mean-P95 historical pr True
sbatch toJobMetricProcess.sh valencia day mean-P95 ssp370 pr True
sbatch toJobMetricProcess.sh valencia day mean-P995 historical pr True
sbatch toJobMetricProcess.sh valencia day mean-P995 ssp370 pr True
sbatch toJobMetricProcess.sh valencia day mean historical pr False
sbatch toJobMetricProcess.sh valencia day mean ssp370 pr False
sbatch toJobMetricProcess.sh valencia day max historical pr False
sbatch toJobMetricProcess.sh valencia day max ssp370 pr False

sbatch toJobMetricProcess.sh valencia 1hr mean-max historical pr True
sbatch toJobMetricProcess.sh valencia 1hr mean-max ssp370 pr True
sbatch toJobMetricProcess.sh valencia 1hr mean-P95 historical pr True
sbatch toJobMetricProcess.sh valencia 1hr mean-P95 ssp370 pr True
sbatch toJobMetricProcess.sh valencia 1hr mean-P995 historical pr True
sbatch toJobMetricProcess.sh valencia 1hr mean-P995 ssp370 pr True
sbatch toJobMetricProcess.sh valencia 1hr mean historical pr False
sbatch toJobMetricProcess.sh valencia 1hr mean ssp370 pr False
sbatch toJobMetricProcess.sh valencia 1hr max historical pr False
sbatch toJobMetricProcess.sh valencia 1hr max ssp370 pr False

sbatch toJobMetricProcess.sh valencia day mean historical tas False
sbatch toJobMetricProcess.sh valencia day mean ssp370 tas False
sbatch toJobMetricProcess.sh valencia 1hr mean historical tas False
sbatch toJobMetricProcess.sh valencia 1hr mean ssp370 tas False

#
sbatch toJobMetricProcess.sh valencia day mean historical pr False
sbatch toJobMetricProcess.sh valencia day mean ssp370 pr False
sbatch toJobMetricProcess.sh valencia 1hr mean historical pr False
sbatch toJobMetricProcess.sh valencia 1hr mean ssp370 pr False
sbatch toJobMetricProcess.sh valencia day max historical pr False
sbatch toJobMetricProcess.sh valencia day max ssp370 pr False
sbatch toJobMetricProcess.sh valencia 1hr max historical pr False
sbatch toJobMetricProcess.sh valencia 1hr max ssp370 pr False

sbatch toJobMetricProcess.sh valencia day max-max historical pr True
sbatch toJobMetricProcess.sh valencia day max-max ssp370 pr True
sbatch toJobMetricProcess.sh valencia 1hr max-max historical pr True
sbatch toJobMetricProcess.sh valencia 1hr max-max ssp370 pr True


sbatch toJobMetricPlot.sh valencia mean-max True True
sbatch toJobMetricPlot.sh valencia max-max True True
sbatch toJobMetricPlot.sh valencia mean-P95 True True
sbatch toJobMetricPlot.sh valencia mean-P995 True True

sbatch toJobMetricPlot.sh valencia mean-max True False
sbatch toJobMetricPlot.sh valencia max-max True False
sbatch toJobMetricPlot.sh valencia mean-P95 True False
sbatch toJobMetricPlot.sh valencia mean-P995 True False



# CMIP5 metrics
sbatch toJobMetricCMIP5.sh rx1day mon-mean Valencia
sbatch toJobMetricCMIP5.sh rx1day yr-mean Valencia
sbatch toJobMetricCMIP5.sh rx1day max-mean Valencia
sbatch toJobMetricCMIP5.sh prhmax mon-mean Valencia
sbatch toJobMetricCMIP5.sh prhmax yr-mean Valencia
sbatch toJobMetricCMIP5.sh prhmax max-mean Valencia

sbatch toJobMetricCMIP5.sh ensemble mon-mean Valencia
sbatch toJobMetricCMIP5.sh ensemble yr-mean Valencia
sbatch toJobMetricCMIP5.sh ensemble max-mean Valencia
 
sbatch toJobMetricCMIP5.sh rx1day mon-mean Iberia
sbatch toJobMetricCMIP5.sh rx1day yr-mean Iberia
sbatch toJobMetricCMIP5.sh rx1day max-mean Iberia
sbatch toJobMetricCMIP5.sh prhmax mon-mean Iberia
sbatch toJobMetricCMIP5.sh prhmax yr-mean Iberia
sbatch toJobMetricCMIP5.sh prhmax max-mean Iberia

sbatch toJobMetricCMIP5.sh ensemble mon-mean Iberia
sbatch toJobMetricCMIP5.sh ensemble yr-mean Iberia
sbatch toJobMetricCMIP5.sh ensemble max-mean Iberia

#PLOTS RELATIVOS fff
sbatch toJobPlotRelCMIP5.sh mon-mean Valencia
sbatch toJobPlotRelCMIP5.sh yr-mean Valencia
sbatch toJobPlotRelCMIP5.sh max-mean Valencia
sbatch toJobPlotRelCMIP5.sh mon-mean Iberia
sbatch toJobPlotRelCMIP5.sh yr-mean Iberia
sbatch toJobPlotRelCMIP5.sh max-mean Iberia
# PLOTS ENSEMBLES
sbatch toJobPlotEnsCMIP5.sh mon-mean Valencia
sbatch toJobPlotEnsCMIP5.sh yr-mean Valencia
sbatch toJobPlotEnsCMIP5.sh max-mean Valencia
sbatch toJobPlotEnsCMIP5.sh mon-mean Iberia
sbatch toJobPlotEnsCMIP5.sh yr-mean Iberia
sbatch toJobPlotEnsCMIP5.sh max-mean Iberia