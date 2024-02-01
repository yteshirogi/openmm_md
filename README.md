# openmm_md
リガンドを含む複合体構造をPDBから取得し、OpenMMによるMDシミュレーションの実行までを自動で行う。
## How to use

### 環境設定

まずファイル類をダウンロードする。
```
$ git clone https://github.com/yteshirogi/openmm_md.git
```
minicondaをインストールし、以下を実行する（Anacondaの場合依存関係の解消で異常に時間がかかったためminiconda推奨）。

```
$ conda create env -n omm -f omm.yml
$ conda activate omm
```
`omm`部分は任意の仮想環境名を設定する。

### 実行
構築したconda仮想環境上で、以下の形式で実行する
```
$ python omm_run.py (pdbid) (ligand_name)
```
例として、3POZで実行する場合は以下のようになる。
```
$ python omm_run.py 3POZ 03P
```
