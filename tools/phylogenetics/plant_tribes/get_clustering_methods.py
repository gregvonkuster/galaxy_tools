import os


def get_clustering_method_options(file_path):
    options = []
    if not os.path.isdir(file_path):
        return options
    methods_dir = os.path.join(file_path, "alns")
    for i, file_name in enumerate(os.listdir(methods_dir)):
        full_path = os.path.join(file_path, file_name)
        options.append((file_name, full_path, i == 0))
    return options
