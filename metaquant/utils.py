def safe_cast_to_list(obj):
    if not isinstance(obj, list):
        return [obj]
    else:
        return obj
