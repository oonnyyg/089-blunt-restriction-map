# blunt-restriction-map

预测线性及环状 DNA 的平端（blunt-end）酶切片段，比较不同酶组合产生的切点与长度分布。仅使用 Python 标准库，不访问任何序列数据库。

## 用法

```
python -m digest_map 输入.json --out 输出目录
```

输出目录中生成：

- `<方案id>.fragments.csv`：每个消化方案一个片段表，列为 `start,length,sequence,wraps_origin`。
- `summary.json`：所有方案的切点（含命中酶与识别起点）、片段数与按长度统计的片段数量直方图。

输入非法时整批失败（退出码 1，错误信息打印到 stderr），且不写任何结果文件。

## 输入格式

```json
{
  "dnas": [
    {"id": "lin", "topology": "linear", "sequence": "AAGATATCGG"},
    {"id": "circ", "topology": "circular", "sequence": "ATCGAT"}
  ],
  "enzymes": [
    {"id": "EcoRV", "recognition": "GATATC"}
  ],
  "digests": [
    {"id": "d1", "dna": "lin", "enzymes": ["EcoRV"]},
    {"id": "d2", "dna": "circ", "enzymes": []}
  ]
}
```

校验规则：

- DNA 序列非空，仅含 ACGT（大小写不敏感，内部统一为大写）；`topology` 为 `linear` 或 `circular`。
- 酶识别序列长度为正偶数，仅含 ACGT，且等于自身反向互补序列；只在中点做平端切割，不模拟甲基化或不完全消化。
- `dnas`、`enzymes`、`digests` 各类的 `id` 分别唯一。
- 每个方案引用一条已定义的 DNA 和无重复的酶集合，空集合允许（表示不切割）。

## 坐标与切割约定

- 采用**零起点碱基间坐标**：位置 `p` 表示第 `p` 个碱基之前的切口，`0` 为 DNA 首碱基之前，长度为 `n` 的线性 DNA 坐标范围是 `0..n`。
- 切点 = 匹配起点 + 识别序列长度的一半；环状 DNA 上对长度 `n` 取模。
- 线性 DNA 不跨端匹配；环状 DNA 允许识别序列跨原点匹配，但识别序列长于 DNA 时该酶无匹配。
- 寻找所有完整匹配，重叠位点不跳过（例如 `ATATAT` 上 `ATAT` 命中起点 0 和 2）。
- 相同位置的切点只切一次，同时保留所有命中酶及其识别起点（见 `summary.json` 的 `hits`），不会因多个酶切在同一位置而产生零长度片段。

## 片段生成

- 线性 DNA：片段包含两端残段，即切点序列 `c1..ck` 产生 `[0,c1), [c1,c2), ..., [ck,n)` 共 `k+1` 段。
- 环状 DNA：从每个切点沿正向前进到下一个切点并跨原点闭合；只有一个切点时产生一条全长线性片段；无切点时保留完整环状分子（单片段，`wraps_origin=false`）。
- 每段报告 `start`（起点坐标）、`length`、`sequence`（正向序列）；仅当 `start + length > n`（序列实际越过原末尾）时 `wraps_origin=true`。
- 所有片段长度之和恒等于原 DNA 长度（程序内部校验）。

## 测试

```
python -m unittest discover -s tests
```
