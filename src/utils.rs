/// TODO
pub fn sync_values(anchor_value: usize, values: &mut [&mut usize]) -> () {
    for value in values {
        **value = anchor_value;
    }
}
