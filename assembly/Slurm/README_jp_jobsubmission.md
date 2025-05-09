# ジョブ投入スクリプトの作成＆実行方法 
for Slurm on NIG Supercomputer, Last updated: 2025-04-19

## ログ保存場所の準備
スパコン上の自分のホームディレクトリに，実行結果の保存場所 `~/log/` ディレクトリを作成しておく。作成は一度だけでよい。
```bash
$ mkdir ~/log/
```
---

## ジョブ投入スクリプトの作成方法
スクリプトは自分のローカルPC上でコードエディタを用いて作成するとよい。（スパコン上で編集したい場合はemacsが使える。）
### スクリプト例
ジョブ実行スクリプト`test.sh`を投入するためのスクリプトの例:
```bash
#!/bin/bash
#SBATCH --output=/home/your_username/log/%x_%j.out.log
#SBATCH --error=/home/your_username/log/%x_%j.err.log
#SBATCH -p epyc
#SBATCH -t 3-23:59:59
#SBATCH -N 1-1 
#SBATCH -n 24
#SBATCH -mem=48G

./test.sh
```

### ヘッダ行（`#`で始まる行）の記述
#### 1. ログ保存フォルダの変更
スクリプトの2~3行目を自分のホームディレクトリ内のログ保存場所のパスに変更（例：`/home/ryohonda/log/%x_%j.out.log`）
#### 2. 実行オプション（リソース要求）を確認・変更
4行目以降で，投入するジョブに必要なスパコン上のリソースを `#SBATCH`オプションで指定して要求できる。実行内容に合わせて要求すること。必要がない場合は指定しなくてもOK。
##### 主なオプション
- `-p` *計算ノード* （例： `epyc`, `short`, `rome`, `medium`）
- `-t` *最大実行時間* (`d-hh:mm:ss`) — 指定しない場合はノードのdefaultのTIMELIMITで自動終了する。(例：*short*ノードは1時間）
- `-N` *最小-最大ノード数* — パラレルジョブ（複数スレッド要求するとき）のみ指定。通常は`-N 1-1`でよい。
- `-n` *CPUスレッド数* — 要求する時は，`-N`も指定すること。
- `--mem=` *ジョブ全体のメモリ量* （例：`--mem=32G`）
- `-J` *ジョブ名（任意）*

### ジョブ実行したいスクリプトの記述
ヘッダ行（#で始まる行）の後に，実行したいスクリプトを記述。スクリプトはパス付きで記述する。実行スクリプトが投入スクリプトと同じディレクトリにある場合は「`./`」を頭につける（`.`=今いるディレクトリを指す）。

次の例のように，複数のスクリプトを記述して順番に実行することも可能。
```bash
#!/bin/bash
#SBATCH --output=/home/your_username/log/%x_%j.out.log
#SBATCH --error=/home/your_username/log/%x_%j.err.log

./test1.sh
./test2.sh
./test3.sh
```
### ジョブ投入スクリプトの保存・転送
拡張子は`.sh`として保存する。改行コードは**必ず「LF」を指定してから保存**する。（Windowsでの改行コードは「CRLF」。この変更忘れるとうまく実行できないので必ず変更すること！！）。

ローカルPCに保存したらFTPでスパコンにアップロード。以降はスパコン上で作業。

---

## パーミッション設定 
### パーミッションの確認
パーミッションの確認は，スクリプトのあるディレクトリに移動して下記シェルコマンドを実行。
```bash
$ ls -l
```
例えば，`-rwxr--r--`　のように最初の4文字目に `x` があれば設定OK（実行権限あり）。
通常のファイルは，`-rw-r--r--`となっていて実行権限(`x`)がない。以下の手順で，作成したスクリプトのパーミッションの変更を行う。

### 実行権限の付加（パーミッションの変更）
スクリプト実行する前にシェルからパーミッション設定をして所有者(`u`)に実行権限(`x`)を付加する。
**ジョブ投入スクリプトとジョブ実行スクリプトの両方**のパーミッションの変更(実行権限の付加）が必要。
```bash
$ chmod u+x スクリプト名.sh
```
まとめて変更したい場合は，スクリプトの置いてあるディレクトリに移動して
```bash
$ chmod u+x *.sh
```
で，`.sh`で終わるファイルすべてに実行権限が付加できる。
終わったら，`ls -l`で実行権限(`x`)が**ジョブ投入スクリプトと実行スクリプトの両方に**正しく変更付加されているかを確認。

## ジョブの投入
*※ジョブ投入は，スパコンのインタラクティブノードから行う（`$ ssh a001` で移動）*

シェルから`sbatch`コマンドで，ジョブ投入スクリプトをジョブスケジューラ(Slurm)に送る。

例：ジョブ投入スクリプト`js_test.sh`の投入
```bash
$ sbatch js_test.sh
```
### ジョブ状況の確認 
投入したジョブの状態を下記コマンドで確認。投入したジョブが`PD`または`R`になっていれば正常に投入できている。
```bash
$ squeue -u 自分のユーザ名
```
ジョブの実行(`R`)が終了すると，リストから消える。ログファイルで実行結果を確認すること（後述）

投入直後なのにジョブがリストにない場合は，投入に失敗した or 実行エラーの可能性がある。ログファイルでエラーを確認すること。

#### ジョブの状態 ####
- **PD** — *PENDING*	実行待ち（リソースが利用可能になるのを待っている）
- **R** — *RUNNING*	実行中
- **CD** — *COMPLETED*	正常終了
- **CA** — *CANCELLED*	ユーザーまたは管理者によってキャンセルされた
- **F** — *FAILED*	エラーで失敗した
- **TO** — *TIMEOUT*	実行時間超過で終了

### ジョブのキャンセル
ジョブのキャンセルは下記シェルコマンドで可能。ジョブ番号は `squeue`（上述）で確認できる。
```bash
$ scancel ジョブ番号
```
### 実行結果・エラーの確認 
実行結果とエラーメッセージは `~/log/`　に保存される。実行後はエラーが出ていないか必ず確認すること。
ログファイルは実行途中でも閲覧可能。進捗状況やエラーが出ていないか適宜確認するとよい。
`log`ディレクトリに移動して`ls -lt`で時間の新しいものから順にログファイルの一覧が表示できる。
```bash
$ cd ~/log/
$ ls -lt
```
#### ログファイルの確認方法
ログファイルは２種類生成される。ログファイルはテキストファイル形式。
- `.out.log` — （出力ログ）スクリプトおよびプログラムによる標準出力（実行状況など）
- `.err.log` — （エラーログ）スパコンやプログラムによるエラーや警告

`ls -lt`の一覧で日付の前の数字がファイルサイズを表す。0以外の場合はなにかしらのログが記録されている。0の場合は空っぽ。
例えば，エラーログのサイズが0の場合，エラーがないことを（ファイルを開かなくても）確認できる。

#### ログファイルの閲覧（moreコマンド）
`more` コマンドを使えば，ローカルPCにダウンロードしなくても内容を確認できる。
```bash
$ more ファイル名
```
ファイルの内容が長い場合は，一番下に白黒反転で 「`--More--(xx%)`」のように表示される。
スペースキーでページ送り。終了する場合は「`q`」をタイプ。

`more`はログファイルだけでなくテキストファイルならなんでも（スクリプトやBLAST出力結果など）閲覧可能なので覚えておくと便利。（`less`という高機能のファイル閲覧コマンドもある）
