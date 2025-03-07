def float_range(value):
    try:
        f = float(value)
        if 0 <= f <= 1:
            return f
        else:
            raise argparse.ArgumentTypeError(f"Value {value} is not between 0 and 1")
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid float value: {value}")