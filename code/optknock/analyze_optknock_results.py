#!/usr/bin/env python3
"""
Analyze OptKnock results across all DAP time points and K values.

Generates:
  - Summary tables
  - Comparison with WT and existing knockout screen
  - Figures (bar charts, heatmaps)
  - LaTeX-formatted results for manuscript

Usage:
    python analyze_optknock_results.py <results_dir>
"""

import sys
import os
import csv
import re
from collections import defaultdict

# Try to import plotting libraries
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("Warning: matplotlib not available. Skipping figure generation.")

try:
    import numpy as np
    HAS_NP = True
except ImportError:
    HAS_NP = False


def read_summary_csv(filepath):
    """Read the merged summary CSV file."""
    results = []
    expected_cols = ['DAP', 'K', 'OptKnock_Lysine', 'Biomass',
                     'Pessimistic_Lysine', 'Optimistic_Lysine',
                     'Biomass_PctWT', 'Growth_Coupled',
                     'ModelStatus', 'SolverStatus', 'Gap']

    with open(filepath, 'r') as f:
        # Skip header (may span multiple lines)
        first_line = f.readline().strip()
        # If header doesn't contain all columns, read next line too
        if first_line.count(',') < 10:
            f.readline()  # skip continuation

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) >= 11:
                row = dict(zip(expected_cols, parts[:11]))
                results.append(row)
    return results


def read_knockouts_csv(filepath):
    """Read the merged knockouts CSV file.
    Handle reaction names containing commas (e.g., R01070[K,c])
    by parsing from both ends of each line."""
    knockouts = []
    with open(filepath, 'r') as f:
        header = f.readline().strip()  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Parse from both ends: DAP,K,...reaction...,WT_Flux,OK_Flux
            parts = line.split(',')
            if len(parts) < 5:
                continue
            dap = parts[0]
            k = parts[1]
            ok_flux = parts[-1]
            wt_flux = parts[-2]
            # Everything in between is the reaction name
            reaction = ','.join(parts[2:-2]).strip().strip('"')
            knockouts.append({
                'DAP': dap,
                'K': k,
                'Reaction': reaction,
                'WT_Flux': wt_flux,
                'OK_Flux': ok_flux
            })
    return knockouts


def print_summary_table(results):
    """Print a formatted summary table."""
    print("\n" + "=" * 110)
    print("  OptKnock Results Summary")
    print("=" * 110)
    print(f"{'DAP':<6} {'K':<4} {'OK_Lysine':>12} {'Biomass':>12} "
          f"{'PessLysine':>12} {'OptiLysine':>12} {'Bio%WT':>8} "
          f"{'Coupled':>8} {'Status':>7}")
    print("-" * 110)

    for row in results:
        dap = row.get('DAP', '?')
        k = row.get('K', '?')
        lys = row.get('OptKnock_Lysine', '0')
        bio = row.get('Biomass', '0')
        pess = row.get('Pessimistic_Lysine', '0')
        opti = row.get('Optimistic_Lysine', '0')
        biopct = row.get('Biomass_PctWT', '?')
        coupled = row.get('Growth_Coupled', '?')
        status = row.get('ModelStatus', '?')

        try:
            lys_f = float(lys)
            bio_f = float(bio)
            pess_f = float(pess)
            opti_f = float(opti)
            biopct_f = float(biopct) if biopct not in ('?', '') else 0
        except ValueError:
            continue

        print(f"{dap:<6} {k:<4} {lys_f:>12.8f} {bio_f:>12.8f} "
              f"{pess_f:>12.8f} {opti_f:>12.8f} {biopct:>7}% "
              f"{coupled:>8} {status:>7}")

    print("=" * 110)


def print_knockout_details(knockouts):
    """Print knockout details organized by DAP and K."""
    print("\n" + "=" * 80)
    print("  Knockout Details")
    print("=" * 80)

    # Group by DAP and K
    grouped = defaultdict(list)
    for row in knockouts:
        key = (row.get('DAP', '?'), row.get('K', '?'))
        grouped[key].append(row)

    for (dap, k), rows in sorted(grouped.items()):
        print(f"\n  DAP={dap}, K={k}:")
        for row in rows:
            rxn = row.get('Reaction', '?').strip()
            wt_flux = row.get('WT_Flux', '0')
            ok_flux = row.get('OK_Flux', '0')
            print(f"    KO: {rxn:<55} WT={wt_flux:>12} OK={ok_flux:>12}")


def find_consistent_knockouts(knockouts, results):
    """Find knockouts that appear consistently across DAP time points."""
    print("\n" + "=" * 80)
    print("  Consistent Knockouts Across DAP Time Points")
    print("=" * 80)

    # Count reaction appearances by K value
    ko_by_k = defaultdict(lambda: defaultdict(set))
    for row in knockouts:
        k = row.get('K', '?')
        rxn = row.get('Reaction', '?').strip()
        dap = row.get('DAP', '?')
        ko_by_k[k][rxn].add(dap)

    for k in sorted(ko_by_k.keys()):
        print(f"\n  K = {k}:")
        # Sort by number of DAPs (most consistent first)
        rxn_counts = [(rxn, daps) for rxn, daps in ko_by_k[k].items()]
        rxn_counts.sort(key=lambda x: -len(x[1]))

        for rxn, daps in rxn_counts:
            dap_str = ', '.join(sorted(daps))
            print(f"    {rxn:<55} [{len(daps)} DAPs: {dap_str}]")


def generate_figures(results, knockouts, output_dir):
    """Generate analysis figures."""
    if not HAS_MPL or not HAS_NP:
        print("Skipping figure generation (matplotlib/numpy not available)")
        return

    fig_dir = os.path.join(output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    # Parse data
    daps = sorted(set(r['DAP'] for r in results))
    ks = sorted(set(int(r['K']) for r in results))

    # ---- Figure 1: Lysine fold-change by DAP and K ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: Lysine flux
    for k in ks:
        k_data = [(r['DAP'], float(r['OptKnock_Lysine']))
                   for r in results if int(r['K']) == k]
        if k_data:
            x_labels = [d[0] for d in k_data]
            y_vals = [d[1] for d in k_data]
            axes[0].plot(x_labels, y_vals, 'o-', label=f'K={k}', markersize=8)

    axes[0].set_xlabel('DAP Time Point', fontsize=12)
    axes[0].set_ylabel('Lysine Flux (R00451[K,p])', fontsize=12)
    axes[0].set_title('OptKnock: Lysine Production', fontsize=14)
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    # Panel B: Biomass % of WT
    for k in ks:
        k_data = [(r['DAP'], float(r['Biomass_PctWT']))
                   for r in results if int(r['K']) == k
                   and r['Biomass_PctWT'] not in ('?', '')]
        if k_data:
            x_labels = [d[0] for d in k_data]
            y_vals = [d[1] for d in k_data]
            axes[1].plot(x_labels, y_vals, 's-', label=f'K={k}', markersize=8)

    axes[1].set_xlabel('DAP Time Point', fontsize=12)
    axes[1].set_ylabel('Biomass (% of WT)', fontsize=12)
    axes[1].set_title('OptKnock: Growth Retention', fontsize=14)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)
    axes[1].axhline(y=100, color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()
    fig.savefig(os.path.join(fig_dir, 'optknock_lysine_biomass.png'),
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fig_dir}/optknock_lysine_biomass.png")

    # ---- Figure 2: Knockout heatmap ----
    if knockouts:
        # Get all unique reactions
        all_rxns = sorted(set(r['Reaction'].strip() for r in knockouts))
        all_daps = sorted(set(r['DAP'] for r in knockouts))
        all_ks = sorted(set(int(r['K']) for r in knockouts))

        # Use the maximum K for the heatmap
        max_k = max(all_ks)
        ko_matrix = np.zeros((len(all_rxns), len(all_daps)))

        for row in knockouts:
            if int(row['K']) == max_k:
                rxn = row['Reaction'].strip()
                dap = row['DAP']
                if rxn in all_rxns and dap in all_daps:
                    ri = all_rxns.index(rxn)
                    di = all_daps.index(dap)
                    ko_matrix[ri, di] = 1

        # Filter to only reactions that appear
        active_rows = np.where(ko_matrix.sum(axis=1) > 0)[0]
        if len(active_rows) > 0:
            fig, ax = plt.subplots(figsize=(8, max(4, len(active_rows) * 0.4)))
            im = ax.imshow(ko_matrix[active_rows], cmap='YlOrRd',
                           aspect='auto', interpolation='nearest')

            ax.set_xticks(range(len(all_daps)))
            ax.set_xticklabels(all_daps, fontsize=10)
            ax.set_yticks(range(len(active_rows)))
            ax.set_yticklabels([all_rxns[i][:40] for i in active_rows],
                               fontsize=8)
            ax.set_xlabel('DAP Time Point', fontsize=12)
            ax.set_title(f'OptKnock Knockouts (K={max_k})', fontsize=14)

            plt.colorbar(im, ax=ax, label='Knocked Out', shrink=0.6)
            plt.tight_layout()
            fig.savefig(os.path.join(fig_dir, 'optknock_knockout_heatmap.png'),
                        dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  Saved: {fig_dir}/optknock_knockout_heatmap.png")

    # ---- Figure 3: Bar chart comparing K values ----
    if len(ks) > 1:
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(daps))
        width = 0.7 / len(ks)

        for idx, k in enumerate(ks):
            k_data = {}
            for r in results:
                if int(r['K']) == k:
                    try:
                        lys = float(r['OptKnock_Lysine'])
                    except (ValueError, KeyError):
                        lys = 0
                    k_data[r['DAP']] = lys

            vals = [k_data.get(d, 0) for d in daps]
            offset = (idx - len(ks)/2 + 0.5) * width
            ax.bar(x + offset, vals, width, label=f'K={k}', alpha=0.85)

        ax.set_xlabel('DAP Time Point', fontsize=12)
        ax.set_ylabel('OptKnock Lysine Flux', fontsize=12)
        ax.set_title('OptKnock: Lysine by K', fontsize=14)
        ax.set_xticks(x)
        ax.set_xticklabels(daps)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        fig.savefig(os.path.join(fig_dir, 'optknock_lysine_foldchange.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fig_dir}/optknock_lysine_foldchange.png")

    print("  Figure generation complete.")


def generate_latex_table(results, knockouts, output_dir):
    """Generate LaTeX table for manuscript."""
    tex_file = os.path.join(output_dir, 'optknock_table.tex')

    with open(tex_file, 'w') as f:
        f.write("% OptKnock Results Table - Auto-generated\n")
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\caption{OptKnock-identified gene knockout strategies for "
                "lysine overproduction in maize endosperm. Lysine flux "
                "(R00451) and biomass retention are shown for different "
                "maximum knockout counts ($K$) across DAP time points.}\n")
        f.write("\\label{tab:optknock}\n")
        f.write("\\small\n")
        f.write("\\begin{tabular}{cc r r r r r}\n")
        f.write("\\hline\n")
        f.write("DAP & $K$ & OptKnock Lysine & Pessimistic Lysine & "
                "Optimistic Lysine & Biomass (\\% WT) & Knocked-Out Reactions \\\\\n")
        f.write("\\hline\n")

        # Group knockouts by (DAP, K)
        ko_by_dk = defaultdict(list)
        for row in knockouts:
            key = (row['DAP'], row['K'])
            ko_by_dk[key].append(row['Reaction'].strip())

        for row in results:
            dap = row.get('DAP', '?')
            k = row.get('K', '?')
            lys = row.get('OptKnock_Lysine', '0')
            pess = row.get('Pessimistic_Lysine', '0')
            opti = row.get('Optimistic_Lysine', '0')
            biopct = row.get('Biomass_PctWT', '?')

            try:
                lys_f = f"{float(lys):.6f}"
                pess_f = f"{float(pess):.6f}"
                opti_f = f"{float(opti):.6f}"
                biopct_f = f"{float(biopct):.1f}" if biopct not in ('?', '') else '?'
            except ValueError:
                lys_f = lys
                pess_f = pess
                opti_f = opti
                biopct_f = biopct

            ko_list = ko_by_dk.get((dap, k), [])
            ko_str = ', '.join([r.replace('_', '\\_') for r in ko_list[:3]])
            if len(ko_list) > 3:
                ko_str += f" (+{len(ko_list)-3} more)"

            f.write(f"{dap} & {k} & {lys_f} & {pess_f} & {opti_f} & "
                    f"{biopct_f}\\% & \\footnotesize{{{ko_str}}} \\\\\n")

        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    print(f"  LaTeX table saved: {tex_file}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_optknock_results.py <results_dir>")
        sys.exit(1)

    results_dir = sys.argv[1]

    summary_file = os.path.join(results_dir, 'all_optknock_summary.csv')
    ko_file = os.path.join(results_dir, 'all_optknock_knockouts.csv')

    if not os.path.exists(summary_file):
        print(f"Error: Summary file not found: {summary_file}")
        sys.exit(1)

    print("=" * 80)
    print("  OptKnock Results Analysis")
    print("=" * 80)

    # Read data
    results = read_summary_csv(summary_file)
    knockouts = read_knockouts_csv(ko_file) if os.path.exists(ko_file) else []

    print(f"  Loaded {len(results)} result rows, {len(knockouts)} knockout entries")

    # Print summary
    print_summary_table(results)

    # Print knockout details
    if knockouts:
        print_knockout_details(knockouts)
        find_consistent_knockouts(knockouts, results)

    # Generate figures
    generate_figures(results, knockouts, results_dir)

    # Generate LaTeX table
    generate_latex_table(results, knockouts, results_dir)

    # Top recommendations
    print("\n" + "=" * 80)
    print("  TOP RECOMMENDATIONS")
    print("=" * 80)

    # Find best results by pessimistic lysine and optimistic lysine
    scored = []
    for row in results:
        try:
            ok_lys = float(row.get('OptKnock_Lysine', '0'))
            pess_lys = float(row.get('Pessimistic_Lysine', '0'))
            opti_lys = float(row.get('Optimistic_Lysine', '0'))
            biopct = float(row.get('Biomass_PctWT', '0'))
            k = int(row.get('K', '0'))
            coupled = row.get('Growth_Coupled', 'no')
            if ok_lys > 0 and biopct > 0:
                scored.append((pess_lys, opti_lys, ok_lys, biopct, k,
                               row['DAP'], coupled))
        except (ValueError, KeyError):
            continue

    # Sort by pessimistic lysine (higher is better = real coupling)
    scored.sort(reverse=True)
    for i, (pess, opti, ok_lys, biopct, k, dap, coupled) in enumerate(scored[:10]):
        print(f"  {i+1}. DAP={dap}, K={k}: "
              f"PessLys={pess:.8f}, OptiLys={opti:.4f}, "
              f"Biomass={biopct:.1f}% WT, Coupled={coupled}")

    print("\n  Analysis complete.")
    print("=" * 80)


if __name__ == '__main__':
    main()
