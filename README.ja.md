# matchfamily — マッチングと重み付け（jamovi モジュール）

交絡がすべて観測されているときの 2 値処置の効果推定：傾向スコアと重なり、共変量バランス、回帰調整、重み付け、二重頑健推定、マッチング、エントロピーバランシング、未観測交絡への感度。 詳しい英語の説明は [README.md](README.md)、中国語の要約は [README.zh.md](README.zh.md) を参照してください。

## 分析（jamovi メニュー：Matching）

| 分析 | 内容 |
|---|---|
| **傾向スコアと重なり** | logit、probit または CBPS。共通サポートとトリミング |
| **共変量バランス** | 標準化平均差、分散比、KS 統計量、Love プロット |
| **重み付けと二重頑健推定量** | RA、IPW、AIPW、IPWRA。ATE/ATT/ATC/ATO |
| **マッチング推定量** | k 近傍、各種距離、完全一致、キャリパー、バイアス調整 |
| **エントロピーバランシング** | 平均（任意で分散・歪度）を厳密に一致 |
| **推定量の比較** | 各手法を並べてフォレスト図 |
| **未観測交絡への感度** | E 値と Rosenbaum バウンド |

各分析の出力には「方法の説明」と「データと設定ガイド」があります。

エンジンは base R で書かれ、Stata（`teffects`、`ebalance`）と R パッケージ（`WeightIt`、`cobalt`、`Matching`、`EValue`、`rbounds`、`PSweight`、`CBPS`）と照合し、推定値と標準誤差は 1e-8 程度で一致します（LaLonde データの IPWRA 標準誤差のみ 1〜2% 差。Stata 自身が未収束を報告）。サンプルデータ：`lalonde`、`sim_confounded.csv`、`sim_overlap.csv`、`sim_hidden.csv`。

## インストール
- **jamovi ライブラリ**：jamovi 右上のモジュールボタンから *jamovi library* を開き、検索してインストール（モジュール公開後）。
- **サイドロード**：`.jmo` ファイルをダウンロードし、同じメニューの *Install from file* から導入します。`.jmo` はビルドした OS、プロセッサ、jamovi のバージョンでのみ動作します。
- **ソースから**：R と `jmvtools` が必要です：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（パスはお使いの jamovi に合わせてください）。

## 引用とライセンス
引用：Sun, M. *matchfamily* (jamovi module). https://github.com/mingchao-sun/matchfamily 。使用した手法の原論文（出力に記載）も併せて引用してください。

ライセンス：`DESCRIPTION` を参照（GPL）。第三者のコードとデータは `NOTICE` に記載。不具合や提案：https://github.com/mingchao-sun/matchfamily/issues
