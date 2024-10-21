def binary_to_decimal(binary_str):
    """Convert a 6-digit binary string to its decimal equivalent."""
    if len(binary_str) != 6 or not all(c in '01' for c in binary_str):
        raise ValueError("Input must be a 6-digit binary string.")
    
    decimal_value = 0
    for i, bit in enumerate(reversed(binary_str)):
        decimal_value += int(bit) * (2 ** i)
    return decimal_value

# Test the function with a sample input
test_binary = "111111"
binary_to_decimal(test_binary)
