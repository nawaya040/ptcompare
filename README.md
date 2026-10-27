# ptcompare — 解釈可能な多変量二標本比較

二つの集団の分布を比較し、全体の差の事後確率と、変数の閾値で説明できる局所的な差を返すRパッケージです。Awaya–Maのhidden Markov Pólya treeの枠組みに基づき、適応的な分割木をSMCで探索し、木を与えた状態推論をmessage passingで計算します。

現在のバージョンは0.1.0です。CRAN提出の準備中です。

## インストール

RcppとC++コンパイラが必要です。WindowsではRtoolsを用意してください。

```r
# install.packages("remotes")
remotes::install_github("nawaya040/ptcompare")
```

## 使い方

```r
library(ptcompare)

set.seed(1)
x <- matrix(rnorm(400), ncol = 2, dimnames = list(NULL, c("marker_A", "marker_B")))
y <- matrix(rnorm(400, mean = 0.8), ncol = 2, dimnames = list(NULL, c("marker_A", "marker_B")))

fit <- compare_pt(x, y)
fit                            # 全体の差の事後確率
summary(fit)                   # 設定、計算の報告、差の大きい分割
regions(fit, prob = 0.9)       # 差のある領域の一覧

plot(fit)                      # 代表木
plot(fit, type = "regions")    # 効果の大きい領域を散布図で並べる
plot(fit, type = "levels")     # 2変数のとき: 木の深さごとに差の確率で塗り分け
```

## Vignette

パッケージの紹介は vignette にまとめています（CRAN提出を見据えて英語）。

```r
vignette("ptcompare")
```

2変量のシミュレーション（仕込んだ差を復元できることを確認）、`mclust::wdbc` の細胞核30変数の解析（サイズでは説明がつかない局所的な差）、差がないデータでの挙動、の3部構成です。

## 自分のデータを比較する

`x` と `y` は、同じ列を同じ順序で持つ数値行列または数値だけのdata.frameです。

```r
library(ptcompare)
set.seed(123)
fit <- compare_pt(x, y)

# 二群のデータをまとめた形式にも対応
fit <- compare_pt(data, group = condition)
```

独立な二標本を対象にしています。NA/Infはエラーとなり、二群を合わせても定数である列は警告付きで除外します。群ラベルがfactorならlevelsの順、その他のベクトルなら最初に現れた順で群1・群2を決めます。

計算予算と統計モデルは、別の引数で指定します。

```r
fit <- compare_pt(
  x, y,
  control = pt_control(particles = 512, effect_draws = 4096),
  prior = pt_prior(max_depth = 8, split_grid = 15, min_node = 20)
)
diagnostics(fit)
```

既定の `pt_prior(midpoint_spike = TRUE)` は論文Appendix Aのspike-and-slab型の分割位置事前です。あるノードを中点で分割すると、その子孫もすべて中点で分割します。すべてのノードで論文の式(2)を使う場合は `midpoint_spike = FALSE` を指定します。

再標本化は論文Section 3.1と同じく、重みのκ乗に比例して祖先を選び、重みをW/aで補正します（既定 `pt_control(resample_power = 0.5)`、`1` で通常の再標本化）。粒子数はSMCの近似精度に関わります。木の深さ、分割候補数、事前確率は推論モデルを変更します。

## 変換の選択

範囲を指定しない `compare_pt(x, y)` は、二群共通の順位変換を既定としています。既知の範囲を `support = cbind(lower, upper)` で指定した場合は、線形変換が既定になります。使用した変換は推定結果に保存され、summaryで表示されます。

```r
fit_rank  <- compare_pt(x, y)                       # 既定: transform = "rank"
fit_paper <- compare_pt(x, y, transform = "range")  # 論文の線形変換
```

`transform = "range"` は論文と同じく、各変数を線形変換で単位区間へ移します。supportを省略すると共通範囲に1%ずつ余裕を加えます。`"rank"` は二群を合わせたmidrankを[0,1]へ移します。群ごとに順位を付ける操作は行いません。変数名、閾値、領域の条件は保存した逆変換によって元の単位で表示されます。同値を保ち、jitterは加えません。

順位変換は元の座標における分割事前と末端分布を変更します。その事後確率は推定した変換に条件付けたもので、変換自体の不確実性を積分していません。

## 結果の読み方

`fit$prob_difference` は有限モデルにおける全体の差の事後確率であり、SMCの木について平均しています。p値ではありません。

`regions()` の各行は代表木の内部ノードです。`condition` で示される領域の中で、`variable <= threshold` の左側へ入る条件付き確率が二群で異なるかを調べています。`prob_difference` はその局所分割の差の事後確率、`effect` は群1対群2の符号付きlog odds ratioの事後平均です。正なら群1の方が左側へ入りやすいという意味になります。`effect_abs` は絶対log odds ratioの事後平均、`lower` と `upper` は差なしの点質量も含む95%信用区間です。

`regions()` は既定で `prob_difference` の大きい順に並びます。効果量の大きい順にする場合は `rank_by = "effect"` を指定します（`plot(fit, type="regions")` と同じ語彙です）。

差の大きい領域の観測そのものを取り出すには `region_rows()` と `region_data()` を使います。

```r
top <- regions(fit, prob = 0.9)$node_id[1]   # 差が最も大きい領域
region_data(fit, top)                        # その領域の観測（群のラベル付き）
mydata[region_rows(fit, top), ]              # 元のdata.frameを同じ条件で絞る
```

`region_rows()` は観測ごとのTRUE/FALSEを、`compare_pt()` に渡した順（`x`の行、続いて`y`の行。`group=`を使った場合は元のdata.frameの行順）で返します。どちらも推定時にデータを保持している必要があります（既定で保持、`pt_control(keep_data = FALSE)` を指定した場合は不可）。

これらの局所量は代表木に条件付けた結果です。代表木はサンプルした候補の中で事後スコアが最大のものを選びます。抽出した領域は親子関係で重なる場合があり、確率閾値による抽出に頻度論的な多重検定保証は付けていません。

図は表示深さで切り詰められ、`+` は下に枝が続くことを示します。図の切り詰めによって推定結果は変わりません。`plot(fit, type="regions")` は差の大きい領域（既定では `prob_difference` が0.5以上で各群50観測以上、最大9個）を並べます。2変数以上では各領域を散布図上の枠として示して枠内の観測を群ごとの記号で描き、1変数では領域内の2群の分布を描きます。`node=` で特定の領域を指定でき、1つだけなら大きな1枚になります。`type="levels"`（2変数のみ）は木の深さごとに領域を差の確率の濃淡で塗ります。

モデル、SMCの目標分布、全体帰無確率、効果量の定義は統計仕様にまとめています。

```r
file.show(system.file("STATISTICAL-SPECIFICATION.md", package = "ptcompare"))
```

## リポジトリから使う場合

以下はソースをクローンした場合の手順で、インストール済みパッケージには含まれません。

**ptcompare.Rproj** を開き、Consoleで次を一度実行するとプロジェクト内の `.local-library` へインストールされます。ソースを修正して再インストールする場合も、RStudioのSession → Restart Rを先に実行してください。

```r
source("scripts/install-local.R")
source("examples/01-local-shift.R")
```

`fit`、`region_table`、`x`、`y` がEnvironmentに残り、Plotsペインに木、局所領域の比較、計算診断が表示されます。Plotsの戻る矢印で前の図を確認できます。

| ファイル | 確認できること |
| --- | --- |
| `examples/00-tryout.R` | 試用スクリプト。基本操作、結果の読み方、iris での実データ例、設定の変え方を順に実行（全体で約10秒） |
| `examples/00-tryout-dimensions.R` | 1変量・2変量・3変量を同じ順番でたどり、次元ごとの見え方を比べる |
| `examples/01-local-shift.R` | 既知の部分集団だけに分布差がある2変量データ |
| `examples/02-dependence-and-null.R` | 周辺分布が同じで相関だけが異なる例（2変量・3変量）と、帰無分布の例 |
| `examples/03-stability.R` | 同じデータ・同じモデルで粒子数と乱数seedを変えて実行 |
| `examples/04-high-dimensional.R` | 20変量の例。既定の順位変換と、論文の線形変換の両方で実行 |

複数回の推定を行う3番目の例は、1回の推定より時間がかかります。実行時間はデータとモデル設定に依存します。

## 検証・再現

これらもリポジトリ内のみです。

- `tests/testthat/`: 独立なR全列挙による固定木・適応木の検証とAPIテスト。
- `evaluation/REPORT.md`: 開発時の検証記録（反復実験、比較手法、計算性能）。
- `evaluation/results/`: seed付きの実験結果CSVと実行環境。

```r
source("scripts/load-local.R")
testthat::test_dir("tests/testthat")

# 評価の再実行。study/rank-studyは既存CSVがあれば未完了分だけ進める。
source("evaluation/run-study.R")
source("evaluation/run-rank-study.R")
source("evaluation/run-verification.R")
source("evaluation/summarize.R")
```

比較実験には追加でMRSパッケージを使います。Windowsでの別プロセスのメモリ測定は `evaluation/run-benchmarks.ps1` にあります。計算核のベンチマークは `evaluation/kernel-benchmark.cpp` です。

参考文献: Awaya, N. and Ma, L. (2024), *Hidden Markov Pólya Trees for High-Dimensional Distributions*, JASA 119(545), 189–201, [doi:10.1080/01621459.2022.2105223](https://doi.org/10.1080/01621459.2022.2105223)。MRSの状態モデルはSoriano and Ma (2017), [doi:10.1111/rssb.12180](https://doi.org/10.1111/rssb.12180) に基づきます。
