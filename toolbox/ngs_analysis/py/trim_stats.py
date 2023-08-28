import pandas as pd
import matplotlib.pyplot as plt
import os

# Snakemake I/O
# INPUTS
list_of_reports = snakemake.input.report
# OUTPUTS
outdir = str(snakemake.output)
# PARAMS
summary_path = str(snakemake.params.summary)
graph_out = str(snakemake.params.graph_summary)

# Create outdir if inexistent
if not os.path.exists(outdir):
	os.makedirs(outdir)

# Imports all trimming reports
for report in list_of_reports:
	sep = os.sep
	sample_name = report.split(sep)[-1]
	col_names = ["category", sample_name]
	df_a = pd.read_csv(report, sep=":\s+", header=None, names=col_names, engine='python')
	try:
		df_increment = pd.read_csv(summary_path, index_col=0)
	except FileNotFoundError:
		df_a.to_csv(summary_path)
		continue
	df_out = pd.concat([df_increment, df_a[sample_name]], join="inner", axis=1)
	df_out.to_csv(summary_path)

# Format dataframes to plot
df_out.index = df_out['category']
df_out = df_out.drop(labels='category', axis=1)
df_out = df_out.T.sort_index()

# Plot trimming summaries
try:
	ax = df_out[['Both Surviving Read Percent',
	             'Forward Only Surviving Read Percent',
	             'Reverse Only Surviving Read Percent',
	             'Dropped Read Percent']].plot.bar(stacked=True).get_figure()
except KeyError:
	pass
	# ax = df_out.plot.bar(stacked=True).get_figure()
try:
	ax.set_size_inches(20, 8)
	ax.savefig(graph_out, dpi=300, bbox_inches='tight')
	plt.close()
except NameError:
	pass
