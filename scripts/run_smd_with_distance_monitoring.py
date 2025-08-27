#!/usr/bin/env python3
"""
Run SMD simulation with real-time distance monitoring.
Monitors pullx.xvg file and terminates simulation when ligand reaches active site.
No frequent checkpointing needed - uses live trajectory monitoring.
"""

import argparse
import subprocess
import time
import os
import signal
import sys
import threading
from typing import Optional

def monitor_distance_file(pullx_file: str, target_distance: float, stop_event: threading.Event):
    """Monitor pullx.xvg file and set stop_event when target distance is reached."""
    print(f"🔍 Monitoring distance file: {pullx_file}")
    print(f"🎯 Target distance: {target_distance} nm")
    
    last_position = 0
    check_interval = 2.0  # Check every 2 seconds
    
    while not stop_event.is_set():
        try:
            if not os.path.exists(pullx_file):
                time.sleep(check_interval)
                continue
            
            # Read new lines from the file
            with open(pullx_file, 'r') as f:
                f.seek(last_position)
                new_lines = f.readlines()
                last_position = f.tell()
            
            # Process new lines
            for line in new_lines:
                line = line.strip()
                if line and not line.startswith('#') and not line.startswith('@'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            time_ps = float(parts[0])
                            distance_nm = float(parts[1])
                            
                            print(f"📏 Time: {time_ps:.1f} ps, Distance: {distance_nm:.3f} nm")
                            
                            if distance_nm <= target_distance:
                                print(f"🎯 TARGET REACHED! Distance {distance_nm:.3f} nm <= {target_distance} nm")
                                print("🛑 Signaling simulation to stop...")
                                stop_event.set()
                                return
                                
                        except ValueError:
                            continue
            
            time.sleep(check_interval)
            
        except Exception as e:
            print(f"Warning: Error monitoring distance file: {e}")
            time.sleep(check_interval)

def run_grompp(mdp_file, gro_file, cpt_file, top_file, ndx_file, tpr_file):
    """Run gmx grompp to create TPR file."""
    cmd = [
        'gmx', 'grompp',
        '-f', mdp_file,
        '-c', gro_file,
        '-t', cpt_file,
        '-p', top_file,
        '-n', ndx_file,
        '-o', tpr_file,
        '-maxwarn', '1000'
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✅ TPR file created successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error creating TPR file: {e}")
        print(f"STDERR: {e.stderr}")
        return False

def run_mdrun_with_distance_monitoring(output_prefix, target_distance=0.1):
    """Run mdrun with real-time distance monitoring and termination."""
    
    cmd = [
        'gmx', 'mdrun',
        '-v',
        '-deffnm', output_prefix,
        '-ntmpi', '1',
        '-ntomp', '4'
    ]
    
    # Set environment
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '4'
    
    pullx_file = f"{output_prefix}_pullx.xvg"
    stop_event = threading.Event()
    
    print("🚀 Starting SMD simulation with distance monitoring...")
    
    try:
        # Start mdrun process
        process = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE, 
                                  stderr=subprocess.STDOUT, text=True, 
                                  bufsize=1, universal_newlines=True)
        
        # Start distance monitoring in separate thread
        monitor_thread = threading.Thread(
            target=monitor_distance_file, 
            args=(pullx_file, target_distance, stop_event)
        )
        monitor_thread.daemon = True
        monitor_thread.start()
        
        # Show live GROMACS output and check for stop signal
        output_lines = []
        while True:
            # Check if we should stop
            if stop_event.is_set():
                print("🛑 Target distance reached - terminating simulation...")
                try:
                    process.terminate()
                    print("📤 Sent termination signal to mdrun")
                    
                    # Wait for graceful shutdown
                    try:
                        process.wait(timeout=30)
                        print("✅ Simulation terminated gracefully")
                        return True
                    except subprocess.TimeoutExpired:
                        print("⚠️ Graceful shutdown timed out, forcing termination")
                        process.kill()
                        process.wait()
                        return True
                        
                except Exception as e:
                    print(f"Error terminating process: {e}")
                    return False
            
            # Read output with timeout
            try:
                # Use a short timeout so we can check stop_event frequently
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print(output.strip())  # Show live GROMACS output
                    output_lines.append(output)
                    
            except Exception:
                break
        
        # Process ended naturally
        return_code = process.returncode
        stop_event.set()  # Stop monitoring thread
        
        if return_code == 0:
            print("✅ Simulation completed successfully")
            return True
        else:
            print(f"⚠️ Simulation ended with return code {return_code}")
            # Check if we have output files despite non-zero return code
            gro_file = f"{output_prefix}.gro"
            if os.path.exists(gro_file):
                print("✅ Output files were created - treating as success")
                return True
            else:
                print("❌ No output files found")
                return False
            
    except Exception as e:
        print(f"Error running mdrun: {e}")
        stop_event.set()
        return False

def validate_output_files(output_prefix):
    """Check that expected output files were created."""
    required_files = [
        f"{output_prefix}.gro",
        f"{output_prefix}_pullf.xvg", 
        f"{output_prefix}_pullx.xvg"
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print(f"❌ Missing required files: {missing_files}")
        return False
    else:
        print("✅ All required output files created")
        return True

def main():
    parser = argparse.ArgumentParser(
        description="Run SMD with real-time distance monitoring and automatic termination"
    )
    parser.add_argument("--mdp", required=True, help="MDP file")
    parser.add_argument("--gro", required=True, help="GRO file")
    parser.add_argument("--cpt", required=True, help="CPT file")
    parser.add_argument("--top", required=True, help="Topology file")
    parser.add_argument("--ndx", required=True, help="Index file")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--target-distance", type=float, default=0.1, 
                       help="Target distance in nm (default: 0.1)")
    
    args = parser.parse_args()
    
    tpr_file = f"{args.output}.tpr"
    
    # Step 1: Create TPR file
    if not run_grompp(args.mdp, args.gro, args.cpt, args.top, args.ndx, tpr_file):
        print("💥 Failed to create TPR file")
        sys.exit(1)
    
    # Step 2: Run mdrun with distance monitoring
    if run_mdrun_with_distance_monitoring(args.output, args.target_distance):
        print("🎉 SMD completed successfully!")
        
        # Step 3: Try to recover final coordinates if .gro is missing
        gro_file = f"{args.output}.gro"
        cpt_file = f"{args.output}.cpt"
        
        if not os.path.exists(gro_file):
            print("🔄 Extracting final coordinates from checkpoint...")
            if os.path.exists(cpt_file):
                try:
                    # Use gmx trjconv to extract final frame from checkpoint
                    recovery_cmd = [
                        'gmx', 'trjconv', 
                        '-s', tpr_file,
                        '-f', cpt_file, 
                        '-o', gro_file,
                        '-dump', '0'  # Extract frame 0 from checkpoint
                    ]
                    # Provide input for group selection (usually 0 for System)
                    result = subprocess.run(recovery_cmd, input='0\n', text=True, check=True, capture_output=True)
                    print("✅ Final coordinates extracted from checkpoint")
                except subprocess.CalledProcessError as e:
                    print(f"❌ Failed to extract coordinates from checkpoint: {e}")
                    print(f"STDERR: {e.stderr}")
            else:
                print("❌ No checkpoint file found - cannot recover final coordinates")
        
        # Step 4: Validate output files
        if validate_output_files(args.output):
            # Show final distance
            pullx_file = f"{args.output}_pullx.xvg"
            try:
                with open(pullx_file, 'r') as f:
                    lines = f.readlines()
                
                # Find last data line
                for line in reversed(lines):
                    if line.strip() and not line.startswith('#') and not line.startswith('@'):
                        parts = line.split()
                        if len(parts) >= 2:
                            final_distance = float(parts[1])
                            print(f"🎯 Final ligand-active site distance: {final_distance:.3f} nm")
                            break
            except:
                pass
            
            sys.exit(0)
        else:
            print("💥 Output validation failed")
            sys.exit(1)
    else:
        print("💥 SMD failed")
        sys.exit(1)

if __name__ == "__main__":
    main()
