# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 23:47:04 2026

@author: ali.s
"""

# ============================================================
# سکشن ۱: ایمپورت‌ها و توابع پایه BPF
# ============================================================
import numpy as np
import pandas as pd
from scipy.integrate import dblquad
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'

def kernel_ts(t, s):
    return t * s

def kernel_t2s(t, s):
    return t**2 * s

def compute_K_matrix(m, kernel):
    h = 1.0 / m
    K = np.zeros((m, m))
    for i in range(m):
        t_low = i * h
        t_high = (i + 1) * h
        for j in range(m):
            s_low = j * h
            s_high = (j + 1) * h
            integral, _ = dblquad(kernel, s_low, s_high,
                                  lambda s: t_low, lambda s: t_high)
            K[i, j] = integral / (h * h)
    return K

def solve_linear(F, K, m):
    h = 1.0 / m
    for i in range(m):
        for j in range(m):
            if K[i, j] < 0:
                F[i, j], F[j, i] = F[j, i], F[i, j]
    I = np.eye(m)
    M = I - h * K
    return np.linalg.solve(M, F)

def reconstruct(Y, t, r, m):
    h = 1.0 / m
    n_idx = min(int(t / h), m-1) if t < 1.0 else m-1
    j_idx = min(int(r / h), m-1) if r < 1.0 else m-1
    return Y[n_idx, j_idx]

# ============================================================
# سکشن ۲: حل‌کننده‌های مثال ۱، ۲ و ۳
# ============================================================

def solve_example1(m):
    f_list = [
        lambda t, r: (1.25 + 0.25 * r) * (-t + np.exp(t)),
        lambda t, r: (1.50 - 0.25 * r) * (-t + np.exp(t)),
        lambda t, r: (1.50 - 0.50 * r) * (-t + np.exp(t)),
        lambda t, r: (1.50 + 0.50 * r) * (-t + np.exp(t))
    ]
    exact_list = [
        lambda t, r: (1.25 + 0.25 * r) * np.exp(t),
        lambda t, r: (1.50 - 0.25 * r) * np.exp(t),
        lambda t, r: (1.50 - 0.50 * r) * np.exp(t),
        lambda t, r: (1.50 + 0.50 * r) * np.exp(t)
    ]
    h = 1.0 / m
    K = compute_K_matrix(m, kernel_ts)
    Y_all = []
    for idx in range(4):
        F = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                F[i, j] = f_list[idx]((i + 0.5) * h, (j + 0.5) * h)
        Y = solve_linear(F, K, m)
        Y_all.append(Y)
    return Y_all, exact_list

def solve_example2(m):
    f_list = [
        lambda t, r: (1.50 + 0.50 * r) * (-5/6 * t**2 + t + 1),
        lambda t, r: (3.00 - 1.00 * r) * (-5/6 * t**2 + t + 1),
        lambda t, r: (2.00 - 1.00 * r) * (-5/6 * t**2 + t + 1),
        lambda t, r: (2.00 + 2.00 * r) * (-5/6 * t**2 + t + 1)
    ]
    exact_list = [
        lambda t, r: (1.50 + 0.50 * r) * (t + 1),
        lambda t, r: (3.00 - 1.00 * r) * (t + 1),
        lambda t, r: (2.00 - 1.00 * r) * (t + 1),
        lambda t, r: (2.00 + 2.00 * r) * (t + 1)
    ]
    h = 1.0 / m
    K = compute_K_matrix(m, kernel_t2s)
    Y_all = []
    for idx in range(4):
        F = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                F[i, j] = f_list[idx]((i + 0.5) * h, (j + 0.5) * h)
        Y = solve_linear(F, K, m)
        Y_all.append(Y)
    return Y_all, exact_list

def g(y):
    return y ** 2

def solve_example3(m):
    a = [0.25, 0.45, 0.35, 0.35]
    b = [0.05, -0.15, -0.10, 0.10]
    integral_const = (np.exp(2) + 1) / 4.0

    def f_source(t, r, idx):
        C = a[idx] + b[idx] * r
        return C * np.exp(t) - t * (C ** 2) * integral_const

    def exact_y(t, r, idx):
        return (a[idx] + b[idx] * r) * np.exp(t)

    h = 1.0 / m
    K = compute_K_matrix(m, kernel_ts)
    Y_all = []

    for idx in range(4):
        F = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                F[i, j] = f_source((i + 0.5) * h, (j + 0.5) * h, idx)

        for i in range(m):
            for j in range(m):
                if K[i, j] < 0:
                    F[i, j], F[j, i] = F[j, i], F[i, j]

        Y = F.copy()
        for _ in range(100):
            Y_new = F + h * (K @ g(Y))
            if np.max(np.abs(Y_new - Y)) < 1e-10:
                Y = Y_new
                break
            Y = Y_new

        Y_all.append(Y)

    exact_list = []
    for idx in range(4):
        exact_list.append(lambda t, r, idx=idx: exact_y(t, r, idx))

    return Y_all, exact_list

# ============================================================
# سکشن ۳: محاسبه خطاها (Full-Grid و Center-Cell)
# ============================================================

def make_table(Y_all, exact_list, t_val, m):
    h = 1.0 / m
    r_points = np.linspace(0, 1, m + 1)
    data = []
    for r in r_points:
        row = [r]
        for idx in range(4):
            approx = reconstruct(Y_all[idx], t_val, r, m)
            exact = exact_list[idx](t_val, r)
            err = abs(approx - exact)
            row.extend([approx, err])
        data.append(row)
    cols = ['r']
    for i in range(4):
        cols.extend([f'App.y{i+1}', f'Err{i+1}'])
    return pd.DataFrame(data, columns=cols)

def compute_max_error(Y_all, exact_list, m):
    h = 1.0 / m
    t_vals = [0.0, 0.125]
    r_vals = np.linspace(0, 1, m + 1)
    max_err = 0
    for t in t_vals:
        for r in r_vals:
            for idx in range(4):
                approx = reconstruct(Y_all[idx], t, r, m)
                exact = exact_list[idx](t, r)
                max_err = max(max_err, abs(approx - exact))
    return max_err

def compute_center_error(Y_all, exact_list, m):
    h = 1.0 / m
    centers = (np.arange(m) + 0.5) * h
    max_errors = np.zeros(4)
    for i, t in enumerate(centers):
        for j, r in enumerate(centers):
            for idx in range(4):
                approx = Y_all[idx][i, j]
                exact = exact_list[idx](t, r)
                err = abs(approx - exact)
                if err > max_errors[idx]:
                    max_errors[idx] = err
    overall_error = np.max(max_errors)
    return max_errors, overall_error

# ============================================================
# سکشن ۴: همگرایی و نرخ‌ها
# ============================================================

def run_example(example_num, m_values):
    print(f"\n{'='*70}")
    print(f"اجرای مثال {example_num}")
    print(f"{'='*70}")
    all_results = {}
    errors = []
    for m in m_values:
        print(f"\n  m = {m} ...")
        if example_num == 1:
            Y_all, exact_list = solve_example1(m)
        elif example_num == 2:
            Y_all, exact_list = solve_example2(m)
        elif example_num == 3:
            Y_all, exact_list = solve_example3(m)
        max_err = compute_max_error(Y_all, exact_list, m)
        errors.append(max_err)
        all_results[m] = (Y_all, exact_list)
        print(f"    Full-Grid Max Error = {max_err:.6e}")
    return all_results, errors

def compute_rate(errors, m_values):
    rates = []
    for i in range(1, len(errors)):
        rate = np.log(errors[i-1] / errors[i]) / np.log(m_values[i] / m_values[i-1])
        rates.append(rate)
    return rates

def compute_center_convergence(example_num, m_values):
    errors = []
    component_errors = []
    for m in m_values:
        if example_num == 1:
            Y_all, exact_list = solve_example1(m)
        elif example_num == 2:
            Y_all, exact_list = solve_example2(m)
        elif example_num == 3:
            Y_all, exact_list = solve_example3(m)
        max_errors, overall_error = compute_center_error(Y_all, exact_list, m)
        errors.append(overall_error)
        component_errors.append(max_errors)
    rates = []
    for i in range(1, len(errors)):
        rate = (np.log(errors[i-1] / errors[i]) /
                np.log(m_values[i] / m_values[i-1]))
        rates.append(rate)
    return errors, component_errors, rates

# ============================================================
# سکشن ۵: جدول مقایسه و کنترل خودکار
# ============================================================

def analyze_all_errors():
    print("\n" + "=" * 70)
    print("مقایسه Full-Grid و Center-Cell Error")
    print("=" * 70)
    m_values = [8, 16, 32, 64]
    summary = {}

    for ex in [1, 2, 3]:
        print(f"\nمثال {ex}:")
        print(f"{'m':<6} {'Full-Grid':<15} {'Center-Cell':<15}")
        print("-" * 40)
        full_errors = []
        center_errors = []
        for m in m_values:
            if ex == 1:
                Y_all, exact_list = solve_example1(m)
            elif ex == 2:
                Y_all, exact_list = solve_example2(m)
            elif ex == 3:
                Y_all, exact_list = solve_example3(m)

            full = compute_max_error(Y_all, exact_list, m)
            _, center = compute_center_error(Y_all, exact_list, m)

            full_errors.append(full)
            center_errors.append(center)

            print(f"{m:<6} {full:<15.6e} {center:<15.6e}")

        summary[ex] = (full_errors, center_errors)

    return summary

def automatic_checks(summary, m_values):
    print("\n" + "=" * 70)
    print("کنترل‌های خودکار")
    print("=" * 70)

    for ex in [1, 2, 3]:
        full_errors, center_errors = summary[ex]
        print(f"\nمثال {ex}:")

        # ۱. آیا خطا با افزایش m کاهش می‌یابد؟
        decreasing = all(full_errors[i] > full_errors[i+1] for i in range(len(full_errors)-1))
        print(f"  ✓ کاهش خطا با m: {'بله ✅' if decreasing else 'خیر ❌'}")

        # ۲. نرخ همگرایی منطقی؟
        rates = []
        for i in range(1, len(full_errors)):
            rate = np.log(full_errors[i-1]/full_errors[i]) / np.log(m_values[i]/m_values[i-1])
            rates.append(rate)
        avg_rate = np.mean(rates)
        last_rate = rates[-1]
        print(f"  ✓ میانگین نرخ: {avg_rate:.2f}")
        print(f"  ✓ آخرین نرخ زوجی: {last_rate:.2f}")

        logical = 0.3 < avg_rate < 3.5
        print(f"  ✓ نرخ منطقی (۰.۳ تا ۳.۵): {'بله ✅' if logical else 'خیر ❌'}")

        # ۳. مقایسه نرخ کلی و آخرین زوج
        diff = abs(avg_rate - last_rate)
        print(f"  ✓ اختلاف نرخ کلی و آخرین زوج: {diff:.2f}")
        print(f"  ✓ پایدار: {'بله ✅' if diff < 0.5 else 'نسبتاً ⚠️'}")
        
# ============================================================
# سکشن اضافه: ذخیره جداول در Excel
# ============================================================
def save_tables_to_excel(example_num, all_results, m_values):
    with pd.ExcelWriter(f'example{example_num}_tables.xlsx') as writer:
        for m in m_values:
            Y_all, exact_list = all_results[m]
            for t_val in [0.0, 0.125]:
                df = make_table(Y_all, exact_list, t_val, m)
                sheet_name = f'm{m}_t{t_val}'
                df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"  ✓ جداول در example{example_num}_tables.xlsx ذخیره شد")        
# ============================================================
# سکشن ۶: رسم نمودارها و اجرای اصلی
# ============================================================

def plot_convergence_all(example_errors, m_values):
    fig, ax = plt.subplots(figsize=(10, 7))

    labels = ['Example 1: Linear, k=ts', 'Example 2: Linear, k=t²s', 'Example 3: Nonlinear, y²']
    colors = ['#e74c3c', '#f39c12', '#2ecc71']
    markers = ['o', 's', 'D']

    for i, (errs, label, color, marker) in enumerate(zip(example_errors, labels, colors, markers)):
        ax.loglog(m_values, errs, marker + '-', color=color, linewidth=2.5,
                  markersize=10, label=label)

    # خطوط مرجع
    ax.loglog(m_values, [0.5/m for m in m_values], '--', color='gray', alpha=0.4, label='O(h)')
    ax.loglog(m_values, [0.5/(m**2) for m in m_values], ':', color='gray', alpha=0.4, label='O(h²)')

    ax.set_xlabel('m (Number of Blocks)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Max Error (log scale)', fontsize=13, fontweight='bold')
    ax.set_title('Convergence Analysis — All Examples (Full-Grid)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig('all_examples_convergence.png', dpi=200)
    plt.savefig('all_examples_convergence.pdf')
    plt.show()
    print("✓ نمودار ترکیبی ذخیره شد")

def plot_center_convergence(center_errors_dict, m_values):
    fig, ax = plt.subplots(figsize=(10, 7))

    labels = ['Example 1', 'Example 2', 'Example 3']
    colors = ['#e74c3c', '#f39c12', '#2ecc71']
    markers = ['o', 's', 'D']

    for ex, (label, color, marker) in zip([1, 2, 3], zip(labels, colors, markers)):
        errs = center_errors_dict[ex]
        ax.loglog(m_values, errs, marker + '-', color=color, linewidth=2.5,
                  markersize=10, label=label)

    ax.set_xlabel('m (Number of Blocks)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Center-Cell Error (log scale)', fontsize=13, fontweight='bold')
    ax.set_title('Center-Cell Convergence — All Examples', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig('center_cell_convergence.png', dpi=200)
    plt.savefig('center_cell_convergence.pdf')
    plt.show()
    print("✓ نمودار Center-Cell ذخیره شد")

def save_report(summary, rates_dict, m_values):
    with open('final_report.txt', 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("گزارش نهایی نتایج\n")
        f.write("=" * 60 + "\n\n")

        f.write("جدول مقایسه بیشینه خطا (Full-Grid):\n")
        f.write(f"{'m':<6} {'Example 1':<15} {'Example 2':<15} {'Example 3':<15}\n")
        f.write("-" * 50 + "\n")
        for i, m in enumerate(m_values):
            f.write(f"{m:<6} {summary[1][0][i]:<15.6e} {summary[2][0][i]:<15.6e} {summary[3][0][i]:<15.6e}\n")

        f.write("\nنرخ‌های همگرایی:\n")
        for ex in [1, 2, 3]:
            rates = rates_dict[ex]
            avg = np.mean(rates)
            f.write(f"  مثال {ex}: نرخ‌ها = {[f'{r:.2f}' for r in rates]}، میانگین = {avg:.2f}\n")

    print("✓ گزارش نهایی در final_report.txt ذخیره شد")

def main():
    print("=" * 70)
    print("حل‌کننده معادلات LIFFIE — سه مثال")
    print("=" * 70)
    print("\nکدام مثال را حل کنم؟")
    print("  ۱: مثال ۱ (خطی، k=ts)")
    print("  ۲: مثال ۲ (خطی، k=t²s)")
    print("  ۳: مثال ۳ (غیرخطی، y²)")
    print("  ۴: هر سه مثال + تحلیل کامل")

    choice = input("\nانتخاب (۱-۴): ").strip()
    m_values = [8, 16, 32, 64]

    if choice in ['1', '2', '3']:
        ex = int(choice)
        results, errors = run_example(ex, m_values)
        rates = compute_rate(errors, m_values)
        print(f"\n  نرخ همگرایی: {np.mean(rates):.2f}")
        save_tables_to_excel(ex, results, [8, 16])
        summary = analyze_all_errors()
        automatic_checks(summary, m_values)
        plot_convergence_all([errors], m_values)
        save_report(summary, {ex: rates}, m_values)

    elif choice == '4':
        all_full_errors = []
        rates_dict = {}
        for ex in [1, 2, 3]:
            results, errors = run_example(ex, m_values)
            all_full_errors.append(errors)
            rates = compute_rate(errors, m_values)
            rates_dict[ex] = rates
            save_tables_to_excel(ex, results, [8, 16])

        summary = analyze_all_errors()
        automatic_checks(summary, m_values)
        plot_convergence_all(all_full_errors, m_values)

        center_errors_dict = {}
        for ex in [1, 2, 3]:
            _, center_errs, _ = compute_center_convergence(ex, m_values)
            center_errors_dict[ex] = center_errs
        plot_center_convergence(center_errors_dict, m_values)

        save_report(summary, rates_dict, m_values)
if __name__ == "__main__":
    main()
    
