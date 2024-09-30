matrix_path='../data/SBS96_De-Novo_Signatures.txt'
output_path='plot'
project='all_virus'
plot_type='96'

python3 -c """
import sigProfilerPlotting as sigPlt
sigPlt.plotSBS(matrix_path='$matrix_path', output_path='$output_path', project='$project', plot_type='$plot_type', percentage=True)
"""