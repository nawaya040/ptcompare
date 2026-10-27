# ptcompare 試用スクリプト: 1変量・2変量・3変量を続けて見る
#
# 目的: CRAN提出に向けて、利用者が既定設定のまま使ったときに何が出てくるかを
# 一度に確かめること。各次元で同じ順番（print → summary → regions → plot）を
# たどるので、次元ごとの見え方の違いを比べられます。
#
# 使い方: ptcompare.Rproj を開き、上から順に Ctrl+Enter で実行してください。
# 右上の Outline（Ctrl+Shift+O）で各セクションへ移動できます。
# Plotsペインの戻る矢印で前の図に戻れます。推定は全部で約1.1秒、
# 通しで実行しても10秒程度です（図は8枚）。
#
# 確認したい点（セクション6にチェックリストがあります）:
#   出力が多すぎないか / 用語が通じるか / 図が説明なしで読めるか。

# 0. 準備 --------------------------------------------------------------------
source("scripts/load-local.R")
packageVersion("ptcompare")
started <- proc.time()[[3]]

# 1. 1変量: 平均も分散も同じで、形だけ違う ------------------------------------
# 群1は標準正規、群2は平均0・分散1の二峰混合。平均の比較では差が見えない例です。
set.seed(101)
n <- 800L
sd_mix <- sqrt(1 - 0.9^2)
x1 <- rnorm(n)
y1 <- ifelse(runif(n) < 0.5, rnorm(n, -0.9, sd_mix), rnorm(n, 0.9, sd_mix))
c(mean_1 = mean(x1), mean_2 = mean(y1), sd_1 = sd(x1), sd_2 = sd(y1))

set.seed(102)
fit1 <- compare_pt(x1, y1)   # 数値ベクトルをそのまま渡せます（列名は X1）
fit1                         # まず全体の差の事後確率
summary(fit1)                # 設定・計算の報告・差の大きい分割
regions(fit1, prob = 0.9, min_n = 30)

plot(fit1)                   # 1変量では同じ変数の閾値が並ぶ木になります
plot(fit1, type = "regions") # 1変量では領域内の2群の分布を重ね描きします
# plot(fit1, type = "levels") は2変量専用なのでエラーになります。

# 2. 2変量: 一部の領域だけに差がある ------------------------------------------
# 群2は marker_A < 0.35 の部分集団だけ marker_B の分布が違います。
set.seed(201)
n <- 1000L
x2 <- cbind(marker_A = runif(n), marker_B = runif(n))
y2 <- cbind(marker_A = runif(n), marker_B = runif(n))
affected <- y2[, "marker_A"] < 0.35
y2[affected, "marker_B"] <- rbeta(sum(affected), 5, 1)

set.seed(202)
fit2 <- compare_pt(x2, y2)
fit2
summary(fit2)
regions(fit2, prob = 0.9, min_n = 40)

plot(fit2)
plot(fit2, type = "regions")  # 2変量では枠がそのまま領域です
plot(fit2, type = "levels")   # 2変量のみ: 深さごとに差の確率で塗り分け

# 3. 3変量: 3次元の箱の中だけ差がある -----------------------------------------
# 群2は X < 0 かつ Y > 0 の部分だけ Z が上にずれています。
set.seed(301)
n <- 900L
draw <- function(n) cbind(X = runif(n, -1, 1), Y = runif(n, -1, 1), Z = rnorm(n))
x3 <- draw(n)
y3 <- draw(n)
in_box <- y3[, "X"] < 0 & y3[, "Y"] > 0
y3[in_box, "Z"] <- y3[in_box, "Z"] + 1.6

set.seed(302)
fit3 <- compare_pt(x3, y3)
fit3
summary(fit3)
regions(fit3, prob = 0.9, min_n = 40)

plot(fit3)
plot(fit3, type = "regions")
# 3変量以上では、枠は描かれている2変数から見た領域です。3つ目の変数で領域の外に
# なる観測は枠の中でも灰色のままなので、枠ではなく記号で所属を読みます。
# 仕込んだ差（Zのずれ）が、別の変数（Y）の条件付き確率の差として出ることがあります。

# 3b. 差の大きい領域の観測を取り出す ------------------------------------------
# regions() は既定で差の事後確率の大きい順。効果量の順にもできます。
head(regions(fit3, prob = 0.5, rank_by = "effect"), 3)

top <- regions(fit3, prob = 0.9)$node_id[1]   # 差が最も大きい領域
sub <- region_data(fit3, top)                 # その領域の観測（群ラベル付き）
dim(sub)
head(sub)
table(sub$group)                              # regions() の n1/n2 と一致します

# region_rows() は観測ごとの TRUE/FALSE。元のdata.frameを絞るのに使えます。
rows <- region_rows(fit3, top)
summary(y3[rows[(nrow(x3) + 1):length(rows)], "Z"])  # 群2の該当部分のZ

# 4. 差がないデータ: 3つの次元で見比べる --------------------------------------
set.seed(401)
null1 <- compare_pt(rnorm(600), rnorm(600))
null2 <- compare_pt(matrix(rnorm(2*600), ncol = 2), matrix(rnorm(2*600), ncol = 2))
null3 <- compare_pt(matrix(rnorm(3*600), ncol = 3), matrix(rnorm(3*600), ncol = 3))
round(c(d1 = null1$prob_difference,
        d2 = null2$prob_difference,
        d3 = null3$prob_difference), 4)
nrow(regions(null2, prob = 0.9))   # 閾値0.9で残る領域の数

# 5. 計算の診断 ---------------------------------------------------------------
plot(fit3, type = "diagnostics")   # ESSの推移（破線=再標本化の閾値、赤点=再標本化）
diagnostics(fit3)[c("final_ess", "resamples", "elapsed_seconds")]

round(sapply(list(d1 = fit1, d2 = fit2, d3 = fit3),
             function(f) diagnostics(f)$elapsed_seconds), 2)
cat(sprintf("このスクリプト全体: %.1f 秒\n", proc.time()[[3]] - started))

# 6. 確認したいこと（CRAN提出前のチェック用メモ）-------------------------------
# 出力について:
#   - print(fit) と summary(fit) の情報量は適切か。summary が長すぎないか。
#   - "Final ESS", "resampling steps" は利用者に意味が伝わるか。
#   - regions() の列（effect, effect_abs, lower, upper, prob_positive）は
#     説明なしで読めるか。幅が広すぎてコンソールで折り返さないか。
# 図について:
#   - 木の図: 既定の max_depth = 3 は適切か。凡例の文言は足りているか。
#   - regions の図: 既定の top = 9 は多すぎないか。1変量の重ね描きは読めるか。
#   - 3変量以上で枠が2変数への射影であることは、図だけで気づけるか。
# エラーメッセージについて:
#   - plot(fit1, type = "levels") の文言は、次に何をすればよいか分かるか。
#   - compare_pt に NA を含むデータを渡したときの文言（下で試せます）。
# 試す場合:
#   bad <- x2; bad[1, 1] <- NA; compare_pt(bad, y2)
