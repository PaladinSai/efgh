import os
import yaml
import importlib.resources
from typing import Any, Dict

class Config:
    """
    用于存储和管理所有配置参数的类，支持属性访问和字典访问。
    Class for storing and managing all configuration parameters, supports attribute and dict access.
    """
    def __init__(self, config_dict: Dict[str, Any]):
        def parse_value(val):
            if isinstance(val, str):
                # 逗号分隔且不是路径（带/或.视为路径/文件名，不分割）
                if "," in val and not any(sep in val for sep in ["/", "."]):
                    return [v.strip() for v in val.split(",") if v.strip()]
                return val.strip()
            elif isinstance(val, dict):
                return Config(val)
            return val

        for k, v in config_dict.items():
            if isinstance(v, dict):
                v = Config(v)
            else:
                v = parse_value(v)
            setattr(self, k, v)
        # 强制 traits, covariates, models 为列表
        if hasattr(self, "gwas"):
            for key in ["traits", "covariates", "models"]:
                val = getattr(self.gwas, key, None)
                if val is None or val == "":
                    setattr(self.gwas, key, [])
                elif isinstance(val, str):
                    setattr(self.gwas, key, [val])
                elif not isinstance(val, list):
                    setattr(self.gwas, key, list(val))
        self._dict = config_dict

    def to_dict(self):
        # 递归转换为普通字典
        # Recursively convert to plain dict
        result = {}
        for k, v in self.__dict__.items():
            if k == "_dict":
                continue
            if isinstance(v, Config):
                result[k] = v.to_dict()
            else:
                result[k] = v
        return result

    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, key, value):
        setattr(self, key, value)

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def load_yaml_from_package(package, resource):
    with importlib.resources.open_text(package, resource, encoding="utf-8") as f:
        return yaml.safe_load(f)

def deep_update(d, u):
    # 递归合并字典
    # Recursively merge dictionaries
    for k, v in u.items():
        if isinstance(v, dict) and k in d and isinstance(d[k], dict):
            deep_update(d[k], v)
        else:
            d[k] = v
    return d

def load_config(default_path=None, user_path=None, cli_args=None, default_pkg=None, default_file=None):
    """
    加载配置，优先级：命令行参数 > 用户yaml > 默认yaml
    Load config, priority: CLI args > user yaml > default yaml
    cli_args: dict，命令行参数（click传入）/ CLI arguments (from click)
    返回Config对象 / Return Config object
    """
    # 1. 加载默认配置 / Load default config
    if default_pkg and default_file:
        config = load_yaml_from_package(default_pkg, default_file)
    elif default_path:
        config = load_yaml(default_path)
    else:
        raise ValueError("No default config source specified")

    # 2. 加载用户自定义配置（如有）/ Load user config if provided
    if user_path and os.path.isfile(user_path):
        user_cfg = load_yaml(user_path)
        config = deep_update(config, user_cfg)

    # 3. 合并命令行参数（如有）/ Merge CLI args if provided
    if cli_args:
        cli_cfg = {k: v for k, v in cli_args.items() if v is not None}
        config = deep_update(config, cli_cfg)

    return Config(config)

def flatten_yaml_dict(d, prefix=""):
    """
    将多层嵌套的yaml配置字典展平成扁平字典，key用.连接
    Flatten nested yaml config dict to flat dict, keys joined by '.'
    用于自动生成命令行参数 / Used for auto-generating CLI options
    """
    items = {}
    for k, v in d.items():
        new_key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            items.update(flatten_yaml_dict(v, new_key))
        else:
            items[new_key] = v
    return items

def get_default_cli_options(default_yaml_path_or_func):
    """
    读取yaml配置，返回所有参数的扁平字典及默认值
    Read yaml config, return flat dict of all parameters and default values
    用于cli.py自动生成命令行参数 / Used for auto-generating CLI options in cli.py
    """
    if callable(default_yaml_path_or_func):
        path = default_yaml_path_or_func()
        with open(path, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
    elif isinstance(default_yaml_path_or_func, str) and os.path.isfile(default_yaml_path_or_func):
        with open(default_yaml_path_or_func, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
    else:
        # 兼容包内资源
        pkg, file = "efgh.configs", "default.yaml"
        with importlib.resources.open_text(pkg, file, encoding="utf-8") as f:
            config = yaml.safe_load(f)
    return flatten_yaml_dict(config)

