import os


def get_clustering_method_options(file_path):
    options = []
    if not os.path.isdir(file_path):
        return options
    methods_dir = os.path.join(file_path, "alns")
    for i, dir_name in enumerate(os.listdir(methods_dir)):
        options.append((dir_name, dir_name, i == 0))
    return options
