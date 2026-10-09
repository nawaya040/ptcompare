# ptcompare — 解釈可能な多変量二標本比較

二つの集団の分布を比較し、全体の差の事後確率と、変数の閾値で説明できる局所的な差を返すRパッケージです。Awaya–Maのhidden Markov Pólya treeの枠組みに基づき、適応的な分割木をSMCで探索し、木を与えた状態推論をmessage passingで計算します。

現在は開発版0.0.1です。CRANへの提出はまだ行っていません。

## RStudioで試す

このフォルダの **ptcompare.Rproj** を開き、Consoleで次を実行してください。このPCではR 4.5.2用に、プロジェクト内の `.local-library` へインストール済みです。

```r
source("examples/01-local-shift.R")
```

`fit`、`region_table`、`x`、`y` がEnvironmentに残り、Plotsペインに木、局所領域の比較、計算診断が表示されます。Plotsの戻る矢印で前の図を確認できます。

```r
summary(fit)
regions(fit, prob = 0.9, min_n = 40)
plot(fit, max_depth = 3)
plot(fit, type = "region", node = 2)
```

別のPCや別のR環境でソースからインストールする場合は、RcppとC++コンパイラが必要です。Windowsでは対応するRtoolsを用意し、Rセッションを再起動してから実行します。

```r
source("scripts/install-local.R")
source("examples/01-local-shift.R")
```

インストール先はプロジェクト内だけです。ソースを修正した後に再インストールする場合も、RStudioのSession → Restart Rを先に実行してください。

## 用意した例

| ファイル | 確認できること |
| --- | --- |
| `examples/00-tryout.R` | 試用スクリプト。基本操作、結果の読み方、iris での実データ例、設定の変え方を順に実行（全体で約10秒） |
| `examples/01-local-shift.R` | 既知の部分集団だけに分布差がある2変量データ。実行はこのPCで約10秒程度 |
| `examples/02-dependence-and-null.R` | 周辺分布が同じで相関だけが異なる例と、帰無分布の例 |
| `examples/03-stability.R` | 同じデータ・同じモデルで粒子数と乱数seedを変えて実行 |
| `examples/04-high-dimensional.R` | 20変量の例。既定の順位変換と、論文の線形変換の両方で実行 |

複数回の推定を行う3番目の例は、1回の推定より時間がかかります。実行時間はデータとモデル設定に依存します。

## 自分のデータを比較する

`x` と `y` は、同じ列を同じ順序で持つ数値行列または数値だけのdata.frameです。

```r
source("scripts/load-local.R")
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

これらの局所量は代表木に条件付けた結果です。代表木はサンプルした候補の中で事後スコアが最大のものを選びます。抽出した領域は親子関係で重なる場合があり、確率閾値による抽出に頻度論的な多重検定保証は付けていません。

図は表示深さで切り詰められ、`+` は下に枝が続くことを示します。図の切り詰めによって推定結果は変わりません。`plot(fit, type="region", node=...)` は、その領域内での群別の分布を元の単位で比較します。

## 検証・再現

- [統計仕様](inst/STATISTICAL-SPECIFICATION.md): モデル、SMCの目標分布、全体帰無確率、効果量の定義。
- [評価レポート](evaluation/REPORT.md): 開発時の検証記録（反復実験、比較手法、計算性能）。
- `tests/testthat/`: 独立なR全列挙による固定木・適応木の検証とAPIテスト。
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
