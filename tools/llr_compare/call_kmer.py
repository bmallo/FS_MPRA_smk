"""Run plasmid-LLR bin/call.py with FiberHMM's hexamer table on a reference containing Ns.

call.py's kmer_arrays() asserts that FiberHMM's encoder marks exactly the reference A/T
positions. On LDLR_MPRA_NotI it does not at 6 positions: the A/T bases flanking the N15
barcode (whose hexamer context contains N) and one base 3 bp from the 3' end. This wrapper
uses the same table-building code but gives those positions zero evidence (no LLR step)
instead of stopping, and prints them. None of them lie in the element +/-100 bp. The
plasmid-LLR repo itself is not modified.

Usage: python tools/llr_compare/call_kmer.py <call.py arguments, with --emissions kmer>
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

CALL = Path('/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/plasmid-llr/bin/call.py')
spec = importlib.util.spec_from_file_location('llr_call', CALL)
mod = importlib.util.module_from_spec(spec)
sys.modules['llr_call'] = mod  # pool workers unpickle _init/_shard by module name
sys.argv = [str(CALL)] + sys.argv[1:]
spec.loader.exec_module(mod)


def kmer_arrays(seq, at_mask, L, model_path=None):
    import fiberhmm
    from fiberhmm.core.bam_reader import encode_from_query_sequence
    from fiberhmm.core.model_io import load_model
    from fiberhmm.inference.tf_recaller import N_CTX, build_llr_tables

    if model_path is None:
        model_path = str(Path(fiberhmm.__file__).parent / 'models/hia5_pacbio.json')
    tab_hit, tab_miss = build_llr_tables(load_model(model_path))
    codes = encode_from_query_sequence(seq, set(range(len(seq))), edge_trim=0, mode='pacbio-fiber', context_size=3)
    info = (codes >= 0) & (codes < N_CTX)
    bad = np.flatnonzero(info != at_mask)
    print(f'[kmer] {len(bad)} A/T positions without a FiberHMM context, given zero evidence (0-based): {bad.tolist()}', flush=True)
    info &= at_mask
    llr_hit, llr_miss = np.zeros(L), np.zeros(L)
    llr_hit[info] = tab_hit[codes[info]]
    llr_miss[info] = tab_miss[codes[info]]
    return llr_hit, llr_miss, info, 0


mod.kmer_arrays = kmer_arrays
mod.main()
