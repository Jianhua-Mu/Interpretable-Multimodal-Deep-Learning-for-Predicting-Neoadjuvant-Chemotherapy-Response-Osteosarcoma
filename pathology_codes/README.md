# Pathomics

这是从 `MIFAPS-main` 中独立出的病理组学工程，只保留全切片图像（WSI）处理、病理特征提取、CLAM/MIL 训练、评估与注意力热图流程。MRI 分类和分割模块、历史检查点、评估结果、热图结果、日志和 Python 缓存均未包含。

## 目录

```text
pathomics/
  data/           数据集与特征读取
  models/         CLAM 和 MIL 模型
  slides/         组织分割、切片与坐标处理
  visualization/  注意力热图
  support/        训练、评估、加载器与序列化
  train.py        病理模型训练入口
  evaluate.py     病理模型评估入口
  generate_splits.py
```

## 数据约定

标注 CSV 至少包含以下列：

- `case_id`：患者或病例标识。
- `slide_id`：切片标识，需对应特征文件名。
- `label`：类别标签；可通过 `--label-column` 指定另一列。

特征根目录应包含 `pt_files/<slide_id>.pt`。标签映射通过 JSON 明确提供，例如 `--label-map '{"negative": 0, "positive": 1}'`。

## 常用命令

在本项目根目录运行：

```powershell
python -m pathomics.generate_splits `
  --dataset-csv D:\data\pathology\labels.csv `
  --output-directory D:\data\pathology\splits `
  --label-map '{"0": 0, "1": 1}'

python -m pathomics.train `
  --dataset-csv D:\data\pathology\labels.csv `
  --feature-root D:\data\pathology\features `
  --split-directory D:\data\pathology\splits `
  --output-directory D:\data\pathology\runs `
  --experiment-name pathology_response `
  --label-map '{"0": 0, "1": 1}'

python -m pathomics.evaluate `
  --dataset-csv D:\data\pathology\labels.csv `
  --feature-root D:\data\pathology\features `
  --model-directory D:\data\pathology\runs\pathology_response_seed82 `
  --output-directory D:\data\pathology\evaluation `
  --label-map '{"0": 0, "1": 1}'
```

运行 `python -m pathomics --help` 可查看全部病理工作流入口。依赖清单见 `requirements.txt`，原始代码许可见 `LICENSE.md`。
