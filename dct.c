#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>
#define PI  3.1415926535897

// gcc -o dct dct.c -lm

void printMatrix(double input[8][8]) {
  int x, y;
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      printf("%8.2f", input[x][y]);
    }
    printf("\n");
  }
}

void printMatrixInt(int input[8][8]) {
  int x, y;
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      printf("%6d", input[x][y]);
    }
    printf("\n");
  }
}

double C(int x) {
  return (x == 0 ? (1 / sqrt(2)) : 1);
}

int main() {
  int x, y, u, v;
  double in[8][8] = {
    {-81,-81,-78,-84,-83,-84,-85,-87},
    {-80,-75,-64,-61,-81,-90,-85,-87},
    {-75,-64,-49,-37,-52,-81,-82,-87},
    {-68,-55,-43,-33,-28,-39,-37,-32},
    {-56,-45,-41,-39,-37,-29,-32,-17},
    {-42,-30,-27,-30,-23,-21,-17,-15},
    {-36,-28,-22,-24,-26,-24,-12,-18},
    {-43,-42,-33,-34,-29,-24,-35,-23}
  };

  int qTable[8][8] = {  
    {16,  11,  10,  16,  24,  40,  51,  61},
    {12,  12,  14,  19,  26,  58,  60,  55},
    {14,  13,  16,  24,  40,  57,  69,  56},
    {14,  17,  22,  29,  51,  87,  80,  62},
    {18,  22,  37,  56,  68, 109, 103,  77},
    {24,  35,  55,  64,  81, 104, 113,  92},
    {49,  64,  78,  87, 103, 121, 120, 101},
    {72,  92,  95,  98, 112, 100, 103,  99}
  };

  double DCT[8][8], sum;
  int outq[8][8];

// DCT MATRIX
  printf("1. Input Matrix: \n");
  printMatrix(in);
  printf("\n\n");

  for (u = 0; u < 8; u++) {
    for (v = 0; v < 8; v++) {

      sum = 0;

      for (x = 0; x < 8; x++) {
        for (y = 0; y < 8; y++) {
          sum += (
            in[x][y] * 
            cos(((2.0 * x + 1) * u * PI) / 16.0) *
            cos(((2.0 * y + 1) * v * PI) / 16.0)
          );
        }
      }

      DCT[u][v] = 0.25 * C(u) * C(v) * sum;
    }
  }

  printf("2. DCT: \n");
  printMatrix(DCT);
  printf("\n");

// QUANTIZED DCT MATRIX
  printf("3. Quantization Table \n");
  printMatrixInt(qTable);
  printf("\n");

  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      outq[x][y] = DCT[x][y] / qTable[x][y];
    }
  }

  printf("4. Quantized DCT Matrix: \n");
  printMatrixInt(outq);
  printf("\n");


// REVERSE QUANTIZATION
  int outReverseQuantization[8][8];

  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      outReverseQuantization[x][y] = outq[x][y] * qTable[x][y];
    }
  }

  printf("5. Reverse Quantization Matrix \n");
  printMatrixInt(outReverseQuantization);
  printf("\n");

// INVERSE DCT MATRIX
  int iDCTMatrix[8][8]; 
  double iSum;

  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++){
      iSum = 0;

      for (u = 0; u < 8; u++) {
        for (v = 0; v < 8; v++) {
          iSum += (
            outReverseQuantization[u][v] * 
            C(u) * 
            C(v) *
            cos((2.0 * x + 1) * u * PI / 16.0) *
            cos((2.0 * y + 1) * v * PI / 16.0)
          );
        }
      }

      iDCTMatrix[x][y] = 0.25 * iSum;
    }
  }

  printf("6. Inverse DCT Matrix \n");
  printMatrixInt(iDCTMatrix);
  printf("\n");

// FIND DIFFERENCE
  int diff[8][8];

  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      diff[x][y] = abs((int) in[x][y] - iDCTMatrix[x][y]);
    }
  }

  printf("7. Difference \n");
  printMatrixInt(diff);
  printf("\n");


// FILE HANDLING
  FILE *fp1;
  fp1 = fopen("matrices.csv", "w");
  if (fp1 == NULL) {
    printf("error");
    return 0;
  }

  // Input Matrix File Write
  fprintf(fp1, "INPUT MATRIX \n");
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      fprintf(fp1, "%.2f%s", in[x][y] + 128, (y < 8 - 1 ? "," : ""));
    }
    fprintf(fp1, "\n");
  }

  fprintf(fp1, "\n");

  fprintf(fp1, "INPUT MATRIX (Reduced) \n");
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      fprintf(fp1, "%.2f%s", in[x][y], (y < 8 - 1 ? "," : ""));
    }
    fprintf(fp1, "\n");
  }

  fprintf(fp1, "\n");

  // Output DCT Matrix File Write
  fprintf(fp1, "DCT MATRIX \n");
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      fprintf(fp1, "%.2f%s", DCT[x][y], (y < 8 - 1 ? "," : ""));
    }
    fprintf(fp1, "\n");
  }

  fprintf(fp1, "\n");

  // Quantization Matrix File Write
  fprintf(fp1, "QUANTIZATION MATRIX \n");
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      fprintf(fp1, "%d%s", qTable[x][y], (y < 8 - 1 ? "," : ""));
    }
    fprintf(fp1, "\n");
  }

  fprintf(fp1, "\n");

  // Quantized DCT Matrix File Write
  fprintf(fp1, "Quantized DCT MATRIX \n");
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      fprintf(fp1, "%d%s", outq[x][y], (y < 8 - 1 ? "," : ""));
    }
    fprintf(fp1, "\n");
  }

  fprintf(fp1, "\n");

  // Quantized DCT Matrix File Write
  fprintf(fp1, "Dequantized DCT MATRIX \n");
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      fprintf(fp1, "%d%s", outReverseQuantization[x][y], (y < 8 - 1 ? "," : ""));
    }
    fprintf(fp1, "\n");
  }

  fprintf(fp1, "\n");

  // Quantized DCT Matrix File Write
  fprintf(fp1, "Inverse the DCT MATRIX to Original \n");
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      fprintf(fp1, "%d%s", iDCTMatrix[x][y], (y < 8 - 1 ? "," : ""));
    }
    fprintf(fp1, "\n");
  }

  fprintf(fp1, "\n");

    // Difference
  fprintf(fp1, "Inverse the DCT MATRIX to Original \n");
  for (x = 0; x < 8; x++) {
    for (y = 0; y < 8; y++) {
      fprintf(fp1, "%d%s", diff[x][y], (y < 8 - 1 ? "," : ""));
    }
    fprintf(fp1, "\n");
  }

  fprintf(fp1, "\n");

  return 0;
}