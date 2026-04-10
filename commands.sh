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