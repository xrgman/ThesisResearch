def normalize_list(data):
    total = sum(data)
    normalized_data = [prob / total for prob in data]

    return normalized_data

