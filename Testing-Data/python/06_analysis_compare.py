import pandas as pd
import numpy as np

df = pd.read_csv(r'd:\Desktop\RTKLIB-VS\Graduation\Testing-Data\data\report\statistics_train_20260510.csv')

print("CSV columns:", list(df.columns))

error_metrics = ['RMS_H', 'RMS_U', 'Mean_H', 'Mean_U', 'Median_H', 'Median_U',
                 'Std_H', 'Std_U', 'P68_H', 'P68_U', 'P95_H', 'P95_U',
                 'RMS.68_H', 'RMS.68_U', 'RMS.95_H', 'RMS.95_U']
ratio_metrics = ['Within_1m_H', 'Within_1m_U', 'Within_3m_H', 'Within_3m_U']
all_metrics = error_metrics + ratio_metrics

# 验证所有列都存在
for m in all_metrics:
    if m not in df.columns:
        print(f"WARNING: column '{m}' not found in CSV")

# 加权平均函数
def w_mean(g, col):
    valid = g[g[col].notna()]
    if len(valid) == 0:
        return np.nan
    return float(np.average(valid[col].values, weights=valid['Count'].values))

# 构建 Dataset x Algorithm 级别数据
rows = []
for (ds, algo), group in df.groupby(['dataset', 'algo']):
    row = {'dataset': ds, 'algo': algo}
    for m in all_metrics:
        row[m] = w_mean(group, m)
    rows.append(row)

df_ds = pd.DataFrame(rows)
print(f"\nDataset x Algorithm 组合数: {len(df_ds)}")
print(f"Dataset数量: {df_ds['dataset'].nunique()}")

# pivot
df_pivot = df_ds.pivot(index='dataset', columns='algo', values=all_metrics)
df_pivot.columns = [f'{m}_{a}' for m, a in df_pivot.columns]
df_pivot = df_pivot.reset_index()

total_ds = len(df_pivot)
print(f"Pivot table行数: {total_ds}\n")

def count_improved(metric, algo_new, algo_base, is_error=True):
    improved = 0
    for _, row in df_pivot.iterrows():
        new_val = row.get(f'{metric}_{algo_new}', np.nan)
        base_val = row.get(f'{metric}_{algo_base}', np.nan)
        if pd.isna(new_val) or pd.isna(base_val):
            continue
        if is_error:
            if new_val < base_val:
                improved += 1
        else:
            if new_val > base_val:
                improved += 1
    return improved

# 保存所有行数据
all_rows = []

for label, algo_new, algo_base in [
    ('SPP-HFUC vs SPP-BRDC', 'SPP-HFUC', 'SPP-BRDC'),
    ('SPP-HFUC-RM vs SPP-BRDC', 'SPP-HFUC-RM', 'SPP-BRDC'),
    ('SPP-HFUC-RM vs SPP-HFUC', 'SPP-HFUC-RM', 'SPP-HFUC'),
]:
    print(f"{'='*65}")
    print(f"【{label}】  (总Dataset数: {total_ds})")
    print(f"{'='*65}")
    print(f"{'指标':<15} | {'提升Dataset数':>8} | {'占比':>8} | {'退化Dataset数':>8}")
    print("-" * 65)

    for m in all_metrics:
        is_err = m in error_metrics
        improved = count_improved(m, algo_new, algo_base, is_error=is_err)
        degraded = total_ds - improved
        pct = improved / total_ds * 100
        print(f"{m:<15} | {improved:>8} / {total_ds} | {pct:>7.1f}% | {degraded:>8}")
        all_rows.append({
            '对比组': label,
            '指标': m,
            '指标类型': '误差类' if is_err else '占比类',
            '提升Dataset数': improved,
            '退化Dataset数': degraded,
            '总Dataset数': total_ds,
            '提升占比(%)': round(pct, 1),
        })

    te = sum(count_improved(m, algo_new, algo_base, is_error=True) for m in error_metrics)
    tr = sum(count_improved(m, algo_new, algo_base, is_error=False) for m in ratio_metrics)
    ta = te + tr
    print("-" * 65)
    print(f"{'误差类小计(16项)':<15} | {te:>8} / {total_ds*16} | {te/total_ds/16*100:>7.1f}% | {total_ds*16-te:>8}")
    print(f"{'占比类小计(4项)':<15} | {tr:>8} / {total_ds*4} | {tr/total_ds/4*100:>7.1f}% | {total_ds*4-tr:>8}")
    print(f"{'综合合计(20项)':<15} | {ta:>8} / {total_ds*20} | {ta/total_ds/20*100:>7.1f}% | {total_ds*20-ta:>8}")
    print()

df_out = pd.DataFrame(all_rows)
df_out.to_csv(r'd:\Desktop\RTKLIB-VS\Graduation\Testing-Data\data\report\compare_results.csv', index=False, encoding='utf-8-sig')
print("结果已保存: compare_results.csv")
