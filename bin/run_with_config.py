#!/usr/bin/env python3
import argparse
import os
import yaml
import subprocess
import shlex
import glob
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Run Kraken ID Parse with a configuration preset')
    parser.add_argument('--preset', required=True, help='Preset configuration to use from kraken_configs.yaml')
    parser.add_argument('--config', help='Path to custom config file', default=None)
    parser.add_argument('--override', action='append', help='Override specific parameters in format key=value', default=[])
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output for debugging')
    return parser.parse_args()

def load_config(config_path=None, verbose=False):
    """Load configuration from specified path or default location"""
    if config_path and os.path.exists(config_path):
        if verbose:
            print(f"Loading config from user-specified path: {config_path}")
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    # Try to find the config file relative to the script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    default_config = os.path.join(repo_root, 'internal', 'kraken_configs.yaml')
    
    if verbose:
        print(f"Script directory: {script_dir}")
        print(f"Repository root: {repo_root}")
        print(f"Looking for default config at: {default_config}")
    
    if os.path.exists(default_config):
        if verbose:
            print(f"Loading config from default path: {default_config}")
        with open(default_config, 'r') as f:
            return yaml.safe_load(f)
    
    raise FileNotFoundError(f"Could not find config file at {default_config}")

def expand_env_vars(value):
    """Expand environment variables in string values"""
    if isinstance(value, str):
        expanded = os.path.expandvars(value)
        return expanded
    return value

def apply_overrides(config, overrides, verbose=False):
    """Apply overrides to the configuration"""
    if verbose:
        print(f"Applying {len(overrides)} overrides: {overrides}")
    
    for override in overrides:
        if '=' not in override:
            print(f"Warning: Invalid override format '{override}', expected key=value")
            continue
            
        key, value = override.split('=', 1)
        orig_value = config.get(key, "(not set)")
        expanded_value = expand_env_vars(value)
        
        # Convert string boolean values to actual booleans
        if value.lower() in ('true', 'false'):
            expanded_value = value.lower() == 'true'
        
        if verbose:
            print(f"Override: {key} = {value}")
            print(f"  Original value: {orig_value}")
            print(f"  Expanded value: {expanded_value}")
        
        config[key] = expanded_value
    
    return config

def main():
    args = parse_args()
    verbose = args.verbose
    
    if verbose:
        print("=== Environment ===")
        print(f"Current working directory: {os.getcwd()}")
        print(f"Python executable: {sys.executable}")
        print(f"Command line arguments: {sys.argv}")
        print("=== Files in current directory ===")
        for f in sorted(os.listdir('.')):
            print(f"  {f}")
    
    try:
        config = load_config(args.config, verbose)
        
        # Get preset configuration
        if args.preset not in config['presets']:
            raise ValueError(f"Preset '{args.preset}' not found in configuration")
        
        if verbose:
            print(f"Using preset: {args.preset}")
            print("Available presets:")
            for preset in config['presets']:
                print(f"  {preset}")
        
        # Start with defaults
        run_config = config.get('defaults', {}).copy()
        
        if verbose:
            print("=== Default configuration ===")
            for key, value in run_config.items():
                print(f"  {key}: {value}")
        
        # Apply preset-specific settings
        run_config.update(config['presets'][args.preset])
        
        if verbose:
            print("=== Configuration after applying preset ===")
            for key, value in run_config.items():
                print(f"  {key}: {value}")
        
        # Apply overrides
        run_config = apply_overrides(run_config, args.override, verbose)
        
        # Expand environment variables in all string values
        if verbose:
            print("=== Expanding environment variables ===")
        
        for key, value in list(run_config.items()):
            expanded = expand_env_vars(value)
            if verbose and isinstance(value, str) and expanded != value:
                print(f"  {key}: {value} -> {expanded}")
            run_config[key] = expanded

        # Drop logo if it doesn't exist to avoid LaTeX failures
        logo_path = run_config.get("logo")
        if isinstance(logo_path, str) and logo_path and not os.path.exists(logo_path):
            if verbose:
                print(f"Logo path not found, removing: {logo_path}")
            run_config.pop("logo", None)
        
        if verbose:
            print("=== Final configuration ===")
            for key, value in run_config.items():
                print(f"  {key}: {value}")
        
        # Build command for kraken_id_parse.py
        script_dir = os.path.dirname(os.path.abspath(__file__))
        kraken_script = os.path.join(script_dir, 'kraken_id_parse.py')
        
        if verbose:
            print(f"Kraken script path: {kraken_script}")
            if not os.path.exists(kraken_script):
                print(f"WARNING: Kraken script not found at {kraken_script}")
        
        # Find FASTQ files using the patterns
        r1_pattern = run_config.get('r1_pattern', '*_R1*.fastq.gz')
        r2_pattern = run_config.get('r2_pattern', '*_R2*.fastq.gz')
        
        if verbose:
            print("=== FASTQ file patterns ===")
            print(f"R1 pattern: {r1_pattern}")
            print(f"R2 pattern: {r2_pattern}")
        
        r1_files = sorted(glob.glob(r1_pattern))
        r2_files = sorted(glob.glob(r2_pattern))
        
        if verbose:
            print("=== Found FASTQ files ===")
            print(f"R1 files ({len(r1_files)}):")
            for f in r1_files:
                print(f"  {f}")
            print(f"R2 files ({len(r2_files)}):")
            for f in r2_files:
                print(f"  {f}")
        
        if not r1_files or not r2_files:
            print(f"WARNING: No FASTQ files found matching patterns:")
            print(f"  R1 pattern: {r1_pattern}")
            print(f"  R2 pattern: {r2_pattern}")
            print(f"Current directory: {os.getcwd()}")
            print("Files in directory:")
            for f in sorted(os.listdir('.')):
                print(f"  {f}")
            print("\nTrying alternative patterns...")
            
            # Try common alternative patterns
            alt_patterns = [
                '*_R1_*.fastq.gz', '*_1.fastq.gz', '*_1_*.fastq.gz', 
                '*_R1.fastq.gz', '*_R1_*.fastq', '*_1.fastq'
            ]
            for pattern in alt_patterns:
                alt_files = glob.glob(pattern)
                if alt_files:
                    print(f"Found files with pattern '{pattern}': {len(alt_files)} files")
                    if not r1_files:  # Only update if we didn't find files with the original pattern
                        r1_files = alt_files
                        run_config['r1_pattern'] = pattern
                        print(f"Using alternative R1 pattern: {pattern}")
            
            alt_patterns = [
                '*_R2_*.fastq.gz', '*_2.fastq.gz', '*_2_*.fastq.gz', 
                '*_R2.fastq.gz', '*_R2_*.fastq', '*_2.fastq'
            ]
            for pattern in alt_patterns:
                alt_files = glob.glob(pattern)
                if alt_files:
                    print(f"Found files with pattern '{pattern}': {len(alt_files)} files")
                    if not r2_files:  # Only update if we didn't find files with the original pattern
                        r2_files = alt_files
                        run_config['r2_pattern'] = pattern
                        print(f"Using alternative R2 pattern: {pattern}")
        
        # Build command parameters
        cmd = [kraken_script]
        
        # Handle special parameters that need different treatment
        special_params = {
            'r1_pattern', 'r2_pattern', 'debug', 'specific'
        }
        
        # Handle file patterns differently - use the actual first file instead of pattern
        if r1_files:
            cmd.extend(['-r1', r1_files[0]])
        elif 'r1_pattern' in run_config:
            cmd.extend(['-r1', run_config['r1_pattern']])
            
        if r2_files:
            cmd.extend(['-r2', r2_files[0]])
        elif 'r2_pattern' in run_config:
            cmd.extend(['-r2', run_config['r2_pattern']])
        
        # Handle debug flag - convert boolean to flag
        if run_config.get('debug'):
            cmd.append('-d')
            if verbose:
                print("Added debug flag (-d)")
        
        # Handle other parameters
        for key, value in run_config.items():
            if key in special_params:
                continue  # Already handled above
                
            # Map config keys to command line arguments
            arg_mapping = {
                'taxon': '-t',
                'kraken_db': '-k', 
                'blast_db': '-b',
                'logo': '-l'
            }
            
            if key in arg_mapping:
                cmd.extend([arg_mapping[key], str(value)])
            else:
                # For any other parameters, use the long form
                cmd.extend([f'--{key}', str(value)])
        
        # Execute the command
        cmd_str = ' '.join(shlex.quote(str(arg)) for arg in cmd)
        print(f"Running command: {cmd_str}")
        
        if verbose:
            print("=== Running subprocess ===")
        
        result = subprocess.run(cmd, capture_output=verbose)
        
        if verbose and result.returncode != 0:
            print("=== Command failed ===")
            print(f"Return code: {result.returncode}")
            print("=== Standard output ===")
            print(result.stdout.decode('utf-8', errors='replace'))
            print("=== Standard error ===")
            print(result.stderr.decode('utf-8', errors='replace'))
        
        return result.returncode
        
    except Exception as e:
        print(f"Error: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == '__main__':
    exit_code = main()
    sys.exit(exit_code)
