#include <stddef.h>
#include "tf_permute.h"


int main() {
    size_t out_array[118];
    simulate_data(234520, .5, 10000000, out_array);
    return 0;
}
