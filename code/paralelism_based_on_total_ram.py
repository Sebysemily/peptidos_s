import os
import psutil

def get_system_ram_limit_mb(leeway_gb=10):
    """
    Get the total virtual memory in MB minus a leeway for the OS and other processes.
    We dynamically scale the leeway (20% of RAM, bounded between 4 GB and leeway_gb GB)
    to ensure smaller systems still have enough memory left to run jobs.
    """
    total_mem_mb = psutil.virtual_memory().total // (1024**2)
    
    # Calculate a dynamic leeway (20% of total RAM, max leeway_gb, min 4GB)
    max_leeway_mb = leeway_gb * 1024
    calculated_leeway_mb = max(4096, min(max_leeway_mb, int(total_mem_mb * 0.20)))
    
    # Remaining memory for Snakemake jobs
    usable_mem_mb = total_mem_mb - calculated_leeway_mb
    return max(2048, usable_mem_mb)  # Ensure at least 2 GB is usable

def set_workflow_resources(workflow, leeway_gb=10):
    """
    Set the global mem_mb resource limit on the Snakemake workflow object.
    """
    usable_mem_mb = get_system_ram_limit_mb(leeway_gb)
    
    if hasattr(workflow, "resource_settings") and hasattr(workflow.resource_settings, "resources"):
        workflow.resource_settings.resources["mem_mb"] = usable_mem_mb
    elif hasattr(workflow, "global_resources"):
        workflow.global_resources["mem_mb"] = usable_mem_mb
    elif hasattr(workflow, "resources") and isinstance(workflow.resources, dict):
        workflow.resources["mem_mb"] = usable_mem_mb

def _fasta_path_from_input(input):
    if hasattr(input, "fasta"):
        return str(input.fasta)
    for inp in input:
        if str(inp).endswith(".fasta"):
            return str(inp)
    if len(input) > 0:
        return str(input[0])
    return None


def _fasta_size_mb(input, fallback_mb=15):
    fasta_path = _fasta_path_from_input(input)
    if fasta_path and os.path.exists(fasta_path):
        try:
            return os.path.getsize(fasta_path) / (1024 * 1024)
        except Exception:
            pass
    return fallback_mb


def get_mem_mb(wildcards, input):
    """
    Estimate memory in MB for a batch job as a function of the input FASTA file size.
    If the file does not exist yet (e.g. during DAG planning/dry-run), fallback to 15 MB.
    """
    file_size_mb = _fasta_size_mb(input, fallback_mb=15)
    # 8.0 GB base memory + 500 MB per MB of FASTA file
    return int(8000 + file_size_mb * 500)


def get_hepad_mem_mb(wildcards, input):
    """
    HEPAD loads pycaret/sklearn pipelines and keeps full feature matrices in RAM.
    Calibrated from OOM kills on 10_25 batches (~2.7 MB FASTA, ~70k seqs):
    peak RSS was ~13-15 GB per job.
    """
    file_size_mb = _fasta_size_mb(input, fallback_mb=2.66)
    # 10 GB model/pipeline overhead + ~2.2 GB per MB of input FASTA
    # → ~16 GB for a typical 10_25 batch (2.7 MB)
    return max(16000, int(10000 + file_size_mb * 2200))
