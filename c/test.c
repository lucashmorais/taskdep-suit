#include <stdio.h>
#include "sha256.h"

int main(int argc, char **argv) {
    SHA256_CTX ctx;
    BYTE buf[SHA256_BLOCK_SIZE];
    BYTE data[] = "Hello!\n";
    sha256_init(&ctx);
    sha256_update(&ctx, data, 7);
    sha256_final(&ctx, buf);
    for(int i = 0; i < SHA256_BLOCK_SIZE; i++) {
        printf("%02x", buf[i]);
    }
    puts("");

    return 0;
}
