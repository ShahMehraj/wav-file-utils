/**
 * @file main.c
 * @brief WAV file processing using FFT
 *
 * This program reads audio samples from a WAV file, performs Fast Fourier Transform (FFT) on the samples,
 * and prints the real and imaginary parts of the FFT result.
 *
 * The program processes the audio samples in blocks of size `BLOCK_SIZE`. Processing in larger block sizes
 * may lead to unresponsiveness, so it is recommended to keep `BLOCK_SIZE` as 1024 or lower.
 *
 * For the remaining samples that are not a multiple of `BLOCK_SIZE`, the program finds the largest power of 4
 * that is less than or equal to the remaining samples and performs FFT on that size.
 *
 *
 * Note: This implementation assumes the availability of the `read_wav` function from the "wave.h" library,
 * as well as the `impeghd_rad2_cplx_fft` function from the "fft.h" library. Make sure to include the correct
 * header files and provide the necessary library files when compiling the program.
 *
 * @date 2023-06-04
 **/

#include <stdio.h>
#include <time.h>
#include "wave.h"
#include "impeghd_type_def.h"
#include "fft.h"

#define BLOCK_SIZE 1024
#include <stdint.h>
#include <string.h>

int main() {
    char file_name[100];
    printf("Enter the name of the WAV file: ");
    scanf("%s", file_name);

    FILE* file = fopen(file_name, "rb");
    if (file == NULL) {
        perror("Failed to open file");
        return 1;
    }

    Header header;
    fread(&header, sizeof(Header), 1, file);

    int option;
    while(1){   
        printf("Choose an option:\n");
        printf("1. File size\n");
        printf("2. Type of format\n");
        printf("3. Number of Channels\n");
        printf("4. Sample Rate\n");
        printf("5. Bits per sample\n");
        printf("6. Size of the data section\n");
        printf("7. FFT of data chunk");
        printf("Enter your choice: ");
        scanf("%d", &option);

        if(option == 7){
            unsigned int num_samples;
            float* ptr_real;
            read_wav(file_name, &num_samples, &ptr_real);
            //num_samples = (num_samples / BLOCK_SIZE) * BLOCK_SIZE;
            float* ptr_imag = calloc(num_samples, sizeof(float));
            float* ptr_scratch = malloc(sizeof(float) * (4 * 1024));

            unsigned int num_of_iterations = num_samples / BLOCK_SIZE;

            for (unsigned int i = 0; i < num_of_iterations; i++) {
                
                unsigned int start_index = i * BLOCK_SIZE;
                unsigned int end_index = (i + 1) * BLOCK_SIZE;
                if(end_index > num_samples) break;

                impeghd_rad2_cplx_fft(
                &ptr_real[start_index], 
                &ptr_imag[start_index], 
                BLOCK_SIZE, 
                ptr_scratch);
            }
            /**
             The following block finds the largest
            power of 4 less than or equal to the
            remaining_samples and computes their fft
            **/
            unsigned int remaining_samples = num_samples % BLOCK_SIZE;
            if (remaining_samples > 0) {
                unsigned int start_index = num_samples - remaining_samples;
                unsigned int end_index = num_samples;

                unsigned int power_of_4 = 1;
                while (power_of_4 * 4 <= remaining_samples) {
                    power_of_4 *= 4;
                }

                impeghd_rad2_cplx_fft(
                    &ptr_real[start_index],
                    &ptr_imag[start_index],
                    power_of_4,
                    ptr_scratch
                );
            }
            for(unsigned int i = 0; i < num_samples; i++){
            printf("%f %c %f i\n", ptr_real[i], (ptr_imag[i] < 0) ? '-' : '+', fabs(ptr_imag[i]));
            }

            free(ptr_real);
            free(ptr_imag);
            free(ptr_scratch);
        }
    

        switch (option) {
            case 1:
                printf("File size: %u bytes\n", header.overall_size);
                break;
            case 2:
                if (header.format_type == 1) {
                    printf("Type of format: PCM\n");
                } else {
                    printf("Type of format: Unknown\n");
                }
                break;
            case 3:
                printf("Number of Channels: %d\n", header.channels);
                break;
            case 4:
                printf("Sample Rate: %u Hz\n", header.sample_rate);
                break;
            case 5:
                printf("Bits per sample: %d bits\n", header.bits_per_sample);
                break;
            case 6:
                printf("Size of the data section: %u bytes\n", header.data_size);
                break;
            default:
                printf("Invalid option\n");
        }
        printf("do you want to continue? [y/n]");
        char response;
        scanf("%c", &response);
        if(response == 'n'){
            break;
        }
        }

    fclose(file);

    return 0;
}


// int main() {
//     const char* file_name = "1khz_Sine_44_1khz.wav";
//     unsigned int num_samples;
//     float* ptr_real;
//     read_wav(file_name, &num_samples, &ptr_real);
//     //num_samples = (num_samples / BLOCK_SIZE) * BLOCK_SIZE;
//     float* ptr_imag = calloc(num_samples, sizeof(float));
//     float* ptr_scratch = malloc(sizeof(float) * (4 * 1024));

//     unsigned int num_of_iterations = num_samples / BLOCK_SIZE;

//     for (unsigned int i = 0; i < num_of_iterations; i++) {
        
//         unsigned int start_index = i * BLOCK_SIZE;
//         unsigned int end_index = (i + 1) * BLOCK_SIZE;
//         if(end_index > num_samples) break;

//         impeghd_rad2_cplx_fft(
//           &ptr_real[start_index], 
//           &ptr_imag[start_index], 
//           BLOCK_SIZE, 
//           ptr_scratch);
//     }
//     /**
//      The following block finds the largest
//      power of 4 less than or equal to the
//      remaining_samples and computes their fft
//     **/
//     unsigned int remaining_samples = num_samples % BLOCK_SIZE;
//     if (remaining_samples > 0) {
//         unsigned int start_index = num_samples - remaining_samples;
//         unsigned int end_index = num_samples;

//         unsigned int power_of_4 = 1;
//         while (power_of_4 * 4 <= remaining_samples) {
//             power_of_4 *= 4;
//         }

//         impeghd_rad2_cplx_fft(
//             &ptr_real[start_index],
//             &ptr_imag[start_index],
//             power_of_4,
//             ptr_scratch
//         );
//     }
//     for(unsigned int i = 0; i < num_samples; i++){
//       printf("%f %c %f i\n", ptr_real[i], (ptr_imag[i] < 0) ? '-' : '+', fabs(ptr_imag[i]));
//     }

//     free(ptr_real);
//     free(ptr_imag);
//     free(ptr_scratch); 

//     return 0;
// }
