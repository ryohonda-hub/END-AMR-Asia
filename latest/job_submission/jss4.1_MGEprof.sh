#!/bin/bash
#SBATCH --output=/home/your_username/log/%x_%j.out
#SBATCH --error=/home/your_username/log/%x_%j.err

./ss4.1_MGEprof.sh

#------- これはジョブ投入用のスクリプトです --------------
# (準備) 実行結果の保存場所 ~/log/ ディレクトリを作成しておく。
#
# 0. ログ保存フォルダの変更：
#		2~3行目を自分のホームディレクトリ内のパスに変更（例：/home/username)
#		*ホームディレクトリ直下に log ディレクトリをあらかじめ作成しておくこと。詳細は下記参照
# 1. ヘッダ行の後を実行したいシェルスクリプト（ジョブ実行スクリプト）に書き換えて保存。
# 2. スクリプトのパーミッション設定を確認（下記参照）
# 3. シェルから下記コマンドでこのジョブ投入用スクリプトを実行
#		$ sbatch このスクリプト名
#
# === ジョブ状況の確認 ===
# 投入したジョブの状況・ジョブ番号は下記シェルコマンドで確認可能。
#		$ squeue -u 自分のユーザ名
# ジョブのキャンセルは下記シェルコマンドで可能。ジョブ番号は squeue で確認。
#		$ scancel ジョブ番号
#
# ** ジョブの状態 **
#	PD	PENDING	実行待ち（リソースが利用可能になるのを待っている）
#	R	RUNNING	実行中
#	CD	COMPLETED	正常終了
#	CA	CANCELLED	ユーザーまたは管理者によってキャンセルされた
#	F	FAILED	エラーで失敗した
#	TO	TIMEOUT	実行時間超過で終了
#
# === 実行結果・エラーの確認 ===
# 実行結果とエラーメッセージは ~/log/　に保存される。実行後はエラーが出ていないか必ず確認すること。
# ホームディレクトリ直下に log ディレクトリをあらかじめ作成しておくこと。
# 		$ mkdir ~/log/
#
# === パーミッション設定 ===
# * スクリプト実行する前にシェルからパーミッション設定をして実行権限(x)を付加してください。一度設定したら次回からは不要。
# ジョブ投入スクリプトとジョブ実行スクリプトの両方のパーミッションの変更(実行権限の付加）が必要です。
# 		$ chmod u+x スクリプト名
# 
# * パーミッションの確認は，スクリプトのあるディレクトリで下記シェルコマンドを実行。
#		$ ls -l
#	例えば，rwxr--r--　のように最初の3文字目に 'x' があれば設定OK（実行権限あり）。
#
# === スクリプトの＃SBATCH ヘッダ行：実行オプション ===
# -n スレッド数
# -t 最大実行時間 (d-hh:mm:ss)。指定しない場合は3日で自動終了する。(shortノードは1時間）
# -J ジョブ名（任意）
