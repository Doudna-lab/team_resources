import pandas as pd
import matplotlib.pyplot as plt
import os


def main():
	# SNAKEMAKE IMPORTS
	# Inputs
	list_of_reports = list(snakemake.input.report)
	# Outputs
	trim_summary_path = str(snakemake.output)
	# Params
	trim_stats_path = str(snakemake.params.summary)
	trim_graph_path = str(snakemake.params.graph_summary)

	if not os.path.exists(trim_summary_path):
		os.makedirs(trim_summary_path)

	for report in list_of_reports:
		sep = os.sep
		sample_name = report.split(sep)[-1]
		col_names = ["category", sample_name]
		df_a = pd.read_csv(report, sep=":\s+", header=None, names=col_names, engine='python')
		try:
			df_increment = pd.read_csv(trim_stats_path, index_col=0)
		except FileNotFoundError:
			df_a.to_csv(trim_stats_path)
			continue

		df_out = pd.concat([df_increment, df_a[sample_name]], join="inner", axis=1)
		df_out.to_csv(trim_stats_path)

	df_out.index = df_out['category']
	df_out = df_out.drop(labels='category', axis=1)
	df_out = df_out.T.sort_index()
	ax = df_out[['Both Surviving Read Percent',
	             'Forward Only Surviving Read Percent',
	             'Reverse Only Surviving Read Percent',
	             'Dropped Read Percent']].plot.bar(stacked=True).get_figure()
	ax.set_size_inches(20, 8)
	ax.savefig(trim_graph_path, dpi=300, bbox_inches='tight')
	plt.close()


if __name__ == "__main__":
	main()
