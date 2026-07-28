import re
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import yt


def extract_pelec_data(filepath):
    """
    Parse a PeleC log file (e.g., DebugNR.log) and extract:
        - time
        - XMOM (x‑momentum)
        - x_velocity MIN and MAX
    from each diagnostic block (each block begins with 'TIME = ... MASS').

    Returns a dictionary with:
        'times'          : list of float
        'xmom'           : list of float
        'xvel_min'       : list of float
        'xvel_max'       : list of float
        'xmom_min'       : float (overall min of XMOM)
        'xmom_max'       : float (overall max)
        'xmom_min_time'  : float (time at which XMOM min occurs)
        'xmom_max_time'  : float (time at which XMOM max occurs)
        'xvel_min_min'   : float (overall min of x_velocity min)
        'xvel_min_max'   : float (overall max of x_velocity min)
        'xvel_max_min'   : float (overall min of x_velocity max)
        'xvel_max_max'   : float (overall max of x_velocity max)
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    blocks = []          # list of dicts: {'time', 'xmom', 'xmin', 'xmax'}
    current_block = None
    in_block = False

    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Start of a new diagnostic block: "TIME = 1.234e-05 MASS = ..."
        if re.match(r'^TIME\s*=\s*[\d.e+-]+\s+MASS', line):
            if in_block and current_block:
                blocks.append(current_block)
            in_block = True
            time_match = re.search(r'TIME\s*=\s*([\d.e+-]+)', line)
            time_val = float(time_match.group(1)) if time_match else None
            current_block = {'time': time_val, 'xmom': None,
                             'xmin': None, 'xmax': None}

        elif in_block and line.startswith('TIME ='):
            # Parse XMOM or x_velocity lines inside the block
            if 'XMOM' in line:
                m = re.search(r'XMOM\s*=\s*([\d.e+-]+)', line)
                if m:
                    current_block['xmom'] = float(m.group(1))
            elif 'x_velocity' in line:
                m_min = re.search(r'MIN\s*=\s*([\d.e+-]+)', line)
                m_max = re.search(r'MAX\s*=\s*([\d.e+-]+)', line)
                if m_min and m_max:
                    current_block['xmin'] = float(m_min.group(1))
                    current_block['xmax'] = float(m_max.group(1))

        else:
            # End of block
            if in_block and current_block:
                # Append only if all fields are present
                if all(v is not None for v in [current_block['time'],
                                                current_block['xmom'],
                                                current_block['xmin'],
                                                current_block['xmax']]):
                    blocks.append(current_block)
                in_block = False
                current_block = None

    # Append last block if any
    if in_block and current_block:
        if all(v is not None for v in [current_block['time'],
                                        current_block['xmom'],
                                        current_block['xmin'],
                                        current_block['xmax']]):
            blocks.append(current_block)

    # Extract into separate lists
    times = [b['time'] for b in blocks]
    xmom = [b['xmom'] for b in blocks]
    xvel_min = [b['xmin'] for b in blocks]
    xvel_max = [b['xmax'] for b in blocks]

    # Compute overall statistics
    def argmin_max(arr):
        if not arr:
            return None, None, None, None
        min_val = min(arr)
        max_val = max(arr)
        min_idx = arr.index(min_val)
        max_idx = arr.index(max_val)
        return min_val, max_val, min_idx, max_idx

    xmom_min, xmom_max, idx_min_xmom, idx_max_xmom = argmin_max(xmom)
    xmin_min, xmin_max, _, _ = argmin_max(xvel_min)
    xmax_min, xmax_max, _, _ = argmin_max(xvel_max)

    result = {
        'times': times,
        'xmom': xmom,
        'xvel_min': xvel_min,
        'xvel_max': xvel_max,
        'xmom_min': xmom_min,
        'xmom_max': xmom_max,
        'xmom_min_time': times[idx_min_xmom] if idx_min_xmom is not None else None,
        'xmom_max_time': times[idx_max_xmom] if idx_max_xmom is not None else None,
        'xvel_min_min': xmin_min,
        'xvel_min_max': xmin_max,
        'xvel_max_min': xmax_min,
        'xvel_max_max': xmax_max,
    }
    return result