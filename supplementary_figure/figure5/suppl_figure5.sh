new_signature="../../data/SBS96_De-Novo_Signatures.txt"
samples="context_96_all_signature_sample.txt"
output="assign_96_denovo_fit"

python3 -c """
from SigProfilerAssignment import Analyzer as Analyze
Analyze.decompose_fit(samples='$samples', signatures='$new_signature',output='$output',input_type='matrix')
"""