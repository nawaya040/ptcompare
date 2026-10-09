# ptcompare 試用スクリプト
#
# 使い方: RStudioで ptcompare.Rproj を開き、このファイルを開いて、
# 上から1行ずつ（または範囲を選択して）Ctrl+Enter で実行してください。
# 右上の「Outline」ボタン（またはCtrl+Shift+O）で各セクションへ移動できます。
# 各セクションは数秒〜十数秒で終わります。

# 0. 準備 ---------------------------------------------------------------
# パッケージのソースを修正した後は、Session > Restart R をしてから
# source("scripts/install-local.R") を一度実行してください。
source("scripts/load-local.R")
packageVersion("ptcompare")

# 1. まず動かす: 一部の領域だけに差がある2変量データ ---------------------
set.seed(1)
n <- 500
x <- data.frame(age = runif(n, 20, 80), score = rnorm(n, 50, 10))
y <- data.frame(age = runif(n, 20, 80), score = rnorm(n, 50, 10))
# 群2では、age < 40 の人だけ score が高い
young <- y$age < 40
y$score[young] <- y$score[young] + 12

set.seed(2)
fit <- compare_pt(x, y)
fit                # 全体の差の事後確率など
summary(fit)       # 局所的な差の上位
plot(fit)          # 代表木（色: 青=群1が左側へ入りやすい、茶=群2が左側へ入りやすい）

# 2. どこに差があるかを調べる -------------------------------------------
reg <- regions(fit, prob = 0.9, min_n = 30)
reg[, c("node_id", "condition", "variable", "threshold", "n1", "n2",
        "prob_difference", "effect", "lower", "upper")]
# condition の領域の中で「variable <= threshold」の側へ入る確率を群間で比べています。
# effect は群1対群2の log odds ratio（正なら群1の方が左側へ入りやすい）。

plot(fit, type = "region")                         # 最も確率の高い領域
plot(fit, type = "region", node = reg$node_id[2])  # 2番目の領域

# 3. 差がないデータでは ---------------------------------------------------
set.seed(3)
x0 <- data.frame(age = runif(n, 20, 80), score = rnorm(n, 50, 10))
y0 <- data.frame(age = runif(n, 20, 80), score = rnorm(n, 50, 10))
set.seed(4)
fit0 <- compare_pt(x0, y0)
fit0$prob_difference
regions(fit0, prob = 0.9)

# 4. 実データの例: iris の2種を4変数で比較 --------------------------------
two <- droplevels(subset(iris, Species != "setosa"))
set.seed(5)
fit_iris <- compare_pt(two[, 1:4], group = two$Species)
summary(fit_iris)
plot(fit_iris)
plot(fit_iris, type = "region")

# 5. 設定を変える ------------------------------------------------------------
# 計算の設定（近似精度）: pt_control()
pt_control()
# 統計モデルの設定: pt_prior()
pt_prior()

# 粒子数を増やす（時間は増えますが、モンテカルロ近似の精度が上がります）
set.seed(6)
fit_big <- compare_pt(x, y, control = pt_control(particles = 1024))
c(default = fit$prob_difference, particles_1024 = fit_big$prob_difference)

# 論文と同じ座標（各変数を線形変換で単位区間へ）を使う
set.seed(7)
fit_range <- compare_pt(x, y, transform = "range")
summary(fit_range)

# 中点スパイク（論文Appendix A）を使わない／通常の再標本化（kappa = 1）にする
set.seed(8)
fit_alt <- compare_pt(x, y, prior = pt_prior(midpoint_spike = FALSE),
                      control = pt_control(resample_power = 1))
fit_alt$prob_difference

# 6. 計算の診断 ------------------------------------------------------------
d <- diagnostics(fit)
d$final_ess          # 最終の有効粒子数
d$elapsed_seconds    # 計算時間（秒）
plot(fit, type = "diagnostics")

# 7. 自分のデータで試す ------------------------------------------------------
# 二群を別々の行列／data.frame で渡す場合（列名と列の順序をそろえる）:
#   fit_my <- compare_pt(data_group1, data_group2)
# 一つの data.frame と群ラベルで渡す場合（ラベルは2種類）:
#   fit_my <- compare_pt(mydata[, c("var1", "var2", "var3")], group = mydata$group)
# 結果の確認:
#   summary(fit_my); regions(fit_my, prob = 0.9); plot(fit_my); plot(fit_my, type = "region")
# 欠測値（NA）を含む行は事前に除いてください。
