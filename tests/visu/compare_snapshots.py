import numpy as np
import visu_ramses
import sys
import os

def compare_snapshots(path1, path2, out1, out2):
    print(f"Comparing Snapshot {out1} from {path1} with Snapshot {out2} from {path2}")

    orig_cwd = os.getcwd()
    
    os.chdir(path1)
    data1 = visu_ramses.load_snapshot(out1)
    
    os.chdir(orig_cwd)
    os.chdir(path2)
    data2 = visu_ramses.load_snapshot(out2)
    
    os.chdir(orig_cwd)

    d1 = data1["data"]
    d2 = data2["data"]

    if len(d1["x"]) != len(d2["x"]):
        print(f"Error: Number of cells mismatch! {len(d1['x'])} vs {len(d2['x'])}")
        return False

    order1 = d1["x"].argsort()
    order2 = d2["x"].argsort()
    # For 2D, we need a better sort. Combine x and y.
    if "y" in d1:
        order1 = np.lexsort((d1["y"], d1["x"]))
        order2 = np.lexsort((d2["y"], d2["x"]))

    vars_to_compare = ["density", "velocity_x", "pressure"]
    if "velocity_y" in d1: vars_to_compare.append("velocity_y")
    
    all_ok = True
    for var in vars_to_compare:
        v1 = d1[var][order1]
        v2 = d2[var][order2]
        
        diff = np.abs(v1 - v2)
        max_diff = np.max(diff)
        avg_diff = np.mean(diff)
        
        print(f"Variable: {var}")
        print(f"  Max Diff: {max_diff:.2e}")
        print(f"  Avg Diff: {avg_diff:.2e}")
        
        if max_diff > 1e-10:
            print(f"  [FAIL] Difference too large!")
            all_ok = False
        else:
            print(f"  [OK] Matches within tolerance.")

    if all_ok:
        print("\nSNAPSHOTS MATCH SUCCESSFULLY!")
    else:
        print("\nSNAPSHOTS DO NOT MATCH!")
    
    return all_ok

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 compare_snapshots.py <dir1> <dir2> <out1> [out2]")
        sys.exit(1)
    
    dir1 = os.path.abspath(sys.argv[1])
    dir2 = os.path.abspath(sys.argv[2])
    out1 = int(sys.argv[3])
    out2 = int(sys.argv[4]) if len(sys.argv) > 4 else out1
    
    sys.path.append(os.path.join(os.path.dirname(__file__), "."))
    compare_snapshots(dir1, dir2, out1, out2)
