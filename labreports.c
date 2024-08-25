// lab -1
/*
#include <stdio.h>

int main()
{
    double arrivalRate = 1.0 / 10 * 60; // customers per hour
    double serviceRate = 1.0 / 5 * 60;  // customers per hour

    // Probability that a customer will not have to wait at the counter
    double P0 = 1 - (arrivalRate / serviceRate);

    // Expected number of customers in the bank
    double ls = arrivalRate / (serviceRate - arrivalRate);

    // Expected time to be spent in the bank
    double ws = 1 / (serviceRate - arrivalRate) * 60;

    printf("Probability that a customer will not have to wait at the counter is %f.\n", P0);
    printf("Expected number of customers in the bank is %f.\n", ls);
    printf("Expected time to be spent in the bank is %f minutes.\n", ws);

    return 0;
}
*/
/*
// lab-2
#include <stdio.h>

int main()
{
    float arrivalRate = 1; // per minute
    float serviceRate = 3; // per minute

    float waitingTime = 1 / (serviceRate - arrivalRate);
    float timeTakenToReachSeat = 1.5; // min

    float totalTimeTaken = waitingTime + timeTakenToReachSeat;

    printf("Average time to get the ticket plus the time to reach the correct seat is %f minutes. \n", totalTimeTaken);
    printf("Hence, the sports fan can expect to be seated for the kick-off.");
    return 0;
}
*/

// lab-3
// 3. WAP to compute PI using Monte Carlo method using C programming language.
/*

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TOTAL_POINTS 1000000

int main()
{
    int inside_circle = 0;
    double x, y;

    // Seed for random number generation
    srand(time(NULL));

    for (int i = 0; i < TOTAL_POINTS; i++)
    {
        // Generate random coordinates in the range [-1, 1]
        x = ((double)rand() / RAND_MAX) * 2 - 1;
        y = ((double)rand() / RAND_MAX) * 2 - 1;

        // Check if the point is inside the unit circle
        if (sqrt(x * x + y * y) <= 1)
            inside_circle++;
    }

    // Approximate π using the ratio of points inside the circle to the total points
    double pi = (double)inside_circle / TOTAL_POINTS * 4;

    printf("Approximated value of PI using Monte Carlo method: %.10f\n", pi);

    return 0;
}
*/

// lab-4
/*
#include <stdio.h>
#include <math.h>

// Function to calculate factorial
unsigned long long factorial(int n)
{
    if (n == 0)
        return 1;
    unsigned long long fact = 1;
    for (int i = 1; i <= n; i++)
    {
        fact *= i;
    }
    return fact;
}

int main()
{
    float arrivalRate = 12.0; // 12 cars/hr
    int x;

    for (x = 0; x <= 15; x++)
    {
        // Poisson probability mass function
        float fx = (pow(exp(1.0), -arrivalRate) * pow(arrivalRate, x)) / factorial(x);
        printf("P(X = %d) = %f\n", x, fx);
    }

    return 0;
}
*/
/*
// lab-5
#include <stdio.h>

int main()
{
    // Transition matrix
    double P[2][2] = {
        {0.4, 0.6}, // P(rain tomorrow | rain today), P(no rain tomorrow | rain today)
        {0.2, 0.8}  // P(rain tomorrow | no rain today), P(no rain tomorrow | no rain today)
    };

    // Matrix to store P * P
    double PP[2][2];

    // Calculate P * P
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            PP[i][j] = 0;
            for (int k = 0; k < 2; k++)
            {
                PP[i][j] += P[i][k] * P[k][j];
            }
        }
    }

    // Initial probabilities
    double initial[2] = {0.0, 1.0}; // Today not raining

    // Probability vector for the day after tomorrow
    double day_after_tomorrow[2];

    // Calculate probability for the day after tomorrow using P * P
    day_after_tomorrow[0] = initial[0] * PP[0][0] + initial[1] * PP[1][0]; // P(rain day after tomorrow)
    day_after_tomorrow[1] = initial[0] * PP[0][1] + initial[1] * PP[1][1]; // P(no rain day after tomorrow)

    // Print the P * P matrix
    printf("Transition matrix P * P:\n");
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            printf("%f ", PP[i][j]);
        }
        printf("\n");
    }

    // Print the result
    printf("Probability that it will not rain the day after tomorrow if it is not raining today: %f\n", day_after_tomorrow[1]);

    return 0;
}
*/

// lab-6
/*
#include <stdio.h>

// Function to print a 2x2 matrix
void printMatrix(double matrix[2][2])
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main()
{
    // Transition matrix P
    double P[2][2] = {
        {0.90, 0.10}, // P(healthy tomorrow | healthy today), P(sick tomorrow | healthy today)
        {0.40, 0.60}  // P(healthy tomorrow | sick today), P(sick tomorrow | sick today)
    };

    // Matrix to store P * P
    double PP[2][2];

    // Calculate P * P
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            PP[i][j] = 0;
            for (int k = 0; k < 2; k++)
            {
                PP[i][j] += P[i][k] * P[k][j];
            }
        }
    }

    // Initial probabilities
    double initial[2] = {0.0, 1.0}; // Sick today

    // Probability vector for the day after tomorrow
    double day_after_tomorrow[2];

    // Calculate probability for the day after tomorrow using P * P
    day_after_tomorrow[0] = initial[0] * PP[0][0] + initial[1] * PP[1][0]; // P(healthy day after tomorrow)
    day_after_tomorrow[1] = initial[0] * PP[0][1] + initial[1] * PP[1][1]; // P(sick day after tomorrow)

    // Print the P * P matrix
    printf("Transition matrix P * P:\n");
    printMatrix(PP);

    // Print the result
    printf("Probability that the person will be sick the day after tomorrow if they are sick today: %.2f\n", day_after_tomorrow[1]);

    return 0;
}
*/
/*
// lab-7
#include <stdio.h>

// Function to print a 2x2 matrix
void printMatrix(double matrix[2][2])
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main()
{
    // Transition matrix P
    double P[2][2] = {
        {0.50, 0.50}, // P(up tomorrow | up today), P(down tomorrow | up today)
        {0.30, 0.70}  // P(up tomorrow | down today), P(down tomorrow | down today)
    };

    // Matrix to store P * P
    double PP[2][2];

    // Calculate P * P
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            PP[i][j] = 0;
            for (int k = 0; k < 2; k++)
            {
                PP[i][j] += P[i][k] * P[k][j];
            }
        }
    }

    // Initial probabilities
    double initial[2] = {0.0, 1.0}; // Market down today

    // Probability vector for the day after tomorrow
    double day_after_tomorrow[2];

    // Calculate probability for the day after tomorrow using P * P
    day_after_tomorrow[0] = initial[0] * PP[0][0] + initial[1] * PP[1][0]; // P(up day after tomorrow)
    day_after_tomorrow[1] = initial[0] * PP[0][1] + initial[1] * PP[1][1]; // P(down day after tomorrow)

    // Print the P * P matrix
    printf("Transition matrix P * P:\n");
    printMatrix(PP);

    // Print the result
    printf("Probability that the market will be down the day after tomorrow if it is down today: %f\n", day_after_tomorrow[1]);

    return 0;
}
*/
/*
// lab-8
#include <stdio.h>

int main()
{
    int m = 100, a = 19, c = 0;
    int X[7];
    X[0] = 63;
    int i;

    printf("The first 7 random numbers are: \n");
    printf("%d\n", X[0]);

    for (i = 0; i < 6; i++)
    {

        X[i + 1] = (a * X[i] + c) % m;

        printf("%d\n", X[i + 1]);
    }
}
*/

/*
// lab-9
#include <stdio.h>

int main()
{
    int m = 100000; // modulus to ensure five-digit integers
    int a = 3;      // multiplier
    int c = 2;      // increment
    int X0 = 4;     // initial seed
    int Xi = X0;    // current value

    printf("Random Integer    Random Variable\n");
    for (int i = 0; i < 10; i++)
    {
        Xi = (a * Xi + c) % m;                                 // Mixed congruential formula
        double random_variable = (double)Xi / m;               // Normalize to get random variable in [0, 1)
        printf("%05d            %.5f\n", Xi, random_variable); // Print results in columns
    }

    return 0;
}
*/
/*

// lab-10
#include <stdio.h>

int main()
{
    int m = 1000; // modulus to ensure three-digit integers
    int a = 3;    // multiplier
    int X0 = 5;   // initial seed
    int Xi = X0;  // current value

    printf("Random Integer   Random Variable\n");
    for (int i = 0; i < 10; i++)
    {
        Xi = (a * Xi) % m;                                       // Multiplicative congruential formula
        double random_variable = (double)Xi / m;                 // Normalize to get random variable in [0, 1)
        printf("%03d              %.3f\n", Xi, random_variable); // Print results in columns
    }

    return 0;
}
*/
// Lab-11
/*
#include <stdio.h>

// Function to compute autocorrelation
void compute_autocorrelation(double *data, int length, double *result)
{
    for (int lag = 0; lag < length; lag++)
    {
        result[lag] = 0.0;
        for (int i = 0; i < length - lag; i++)
        {
            result[lag] += data[i] * data[i + lag];
        }
        // Normalize by the number of points
        result[lag] /= (length - lag);
    }
}

int main()
{
    // Example data
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    int length = sizeof(data) / sizeof(data[0]);

    // Array to hold autocorrelation results
    double result[length];

    // Compute autocorrelation
    compute_autocorrelation(data, length, result);

    // Print results
    printf("Autocorrelation:\n");
    for (int lag = 0; lag < length; lag++)
    {
        printf("Lag %d: %.2f\n", lag, result[lag]);
    }

    return 0;
}
*/
/*
// lab-12

#include <stdio.h>

#include <math.h>
#include <stdlib.h>

// Function to calculate the theoretical CDF
double theoretical_cdf(double x)
{
    return 1 - pow(0.9, x + 1);
}

// Function to sort an array
void sort_array(int *array, int size)
{
    for (int i = 0; i < size - 1; i++)
    {
        for (int j = i + 1; j < size; j++)
        {
            if (array[i] > array[j])
            {
                int temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }
    }
}

// Function to calculate the empirical CDF
double empirical_cdf(int *gaps, int size, int x)
{
    int count = 0;
    for (int i = 0; i < size; i++)
    {
        if (gaps[i] <= x)
        {
            count++;
        }
    }
    return (double)count / size;
}

int main()
{
    int n, target, max_gap;
    double alpha, D_alpha;

    // Input sequence size
    printf("Enter the number of elements in the sequence: ");
    scanf("%d", &n);

    int sequence[n];

    // Input sequence values
    printf("Enter the sequence values: ");
    for (int i = 0; i < n; i++)
    {
        scanf("%d", &sequence[i]);
    }

    // Input target value
    printf("\nEnter the target value: ");
    scanf("%d", &target);

    // Input maximum gap to count
    printf("Enter the maximum gap to count: ");
    scanf("%d", &max_gap);

    // Input significance level alpha
    printf("Enter the significance level alpha (e.g., 0.05): ");
    scanf("%lf", &alpha);

    // Input critical value D_alpha from KS table
    printf("Enter the critical value D_alpha: ");
    scanf("%lf", &D_alpha);

    int gaps[n];
    int gap_counts[max_gap + 1];

    // Initialize gaps and gap_counts arrays to 0
    for (int i = 0; i < n; i++)
    {
        gaps[i] = 0;
    }
    for (int i = 0; i <= max_gap; i++)
    {
        gap_counts[i] = 0;
    }

    // Perform gap test
    int gap = 0;
    int gap_index = 0;
    for (int i = 0; i < n; i++)
    {
        if (sequence[i] == target)
        {
            if (gap > 0)
            {
                if (gap <= max_gap)
                {
                    gap_counts[gap]++;
                }
                else
                {
                    gap_counts[max_gap]++;
                }
                gaps[gap_index++] = gap;
                gap = 0;
            }
        }
        else
        {
            gap++;
        }
    }

    // Sort the gaps array
    sort_array(gaps, gap_index);

    // Calculate D
    double D = 0.0;
    printf("\nCalculations:\n");
    printf(" x |   F(x)   |  S_N(x)  | |F(x) - S_N(x)|\n");
    printf("---|----------|----------|---------------\n");
    for (int i = 0; i <= max_gap; i++)
    {
        double F_x = theoretical_cdf(i);
        double S_N_x = empirical_cdf(gaps, gap_index, i);
        double D_current = fabs(F_x - S_N_x);
        printf("%2d   %8.6lf   %8.6lf   %13.6lf\n", i, F_x, S_N_x, D_current);
        if (D_current > D)
        {
            D = D_current;
        }
    }

    // Output the maximum deviation D
    printf("\nMaximum Deviation D: %lf\n", D);

    // Determine whether to reject the null hypothesis
    if (D < D_alpha)
    {
        printf("D_cal < D_alpha, Null hypothesis is not rejected.\n");
    }
    else
    {
        printf("D_cal >= D_alpha, Null hypothesis is rejected.\n");
    }

    return 0;
}
*/
/*
// lab-14

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to calculate the Chi-Square statistic
double chi_square_statistic(int *observed, double *expected, int size)
{
    double chi_square = 0.0;
    for (int i = 0; i < size; i++)
    {
        double difference = observed[i] - expected[i];
        chi_square += (difference * difference) / expected[i];
    }
    return chi_square;
}

int main()
{
    int num_categories;
    double alpha, critical_value;

    // Input number of categories
    printf("Enter the number of categories: ");
    scanf("%d", &num_categories);

    int observed[num_categories];
    double expected[num_categories];

    // Input observed frequencies
    printf("Enter the observed frequencies:\n");
    for (int i = 0; i < num_categories; i++)
    {
        scanf("%d", &observed[i]);
    }

    // Input expected frequencies
    printf("Enter the expected frequencies:\n");
    for (int i = 0; i < num_categories; i++)
    {
        scanf("%lf", &expected[i]);
    }

    // Input significance level alpha
    printf("Enter the significance level alpha (e.g., 0.05): ");
    scanf("%lf", &alpha);

    // Input critical value from Chi-Square distribution table
    printf("Enter the critical value from the Chi-Square distribution table: ");
    scanf("%lf", &critical_value);

    // Calculate Chi-Square statistic
    double chi_square = chi_square_statistic(observed, expected, num_categories);

    // Output calculations
    printf("\nCalculations:\n");
    printf(" Category | Observed | Expected | (Observed - Expected)^2 / Expected\n");
    printf("----------|----------|----------|-----------------------------------\n");

    for (int i = 0; i < num_categories; i++)
    {
        double difference = observed[i] - expected[i];
        double chi_square_contribution = (difference * difference) / expected[i];
        printf("%9d | %8d | %8.2lf | %30.2lf\n", i + 1, observed[i], expected[i], chi_square_contribution);
    }

    printf("\nChi-Square Statistic: %lf\n", chi_square);

    // Determine whether to reject the null hypothesis
    if (chi_square < critical_value)
    {
        printf("Chi-Square statistic < Critical value, Null hypothesis is not rejected.\n");
    }
    else
    {
        printf("Chi-Square statistic >= Critical value, Null hypothesis is rejected.\n");
    }

    return 0;
}
* /
    /*
    // lab-13
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>

    // Function to compare two doubles for qsort
    int compare(const void *a, const void *b)
    {
        return (*(double *)a > *(double *)b) - (*(double *)a < *(double *)b);
    }

    // Function to perform the KS test
    double ks_test(double *sample, int n)
    {
        qsort(sample, n, sizeof(double), compare);

        double d_plus = 0.0, d_minus = 0.0;

        printf("| %4s | %6s | %8s | %8s | %15s | %15s |\n", "i", "Ri", "i/n", "(i-1)/n", "D+ = i/n - Ri", "D- = Ri - (i-1)/n");
        printf("|------|--------|----------|----------|-----------------|-----------------|\n");

        for (int i = 0; i < n; i++)
        {
            double ri = sample[i];
            double d_plus_i = (i + 1.0) / n - ri;
            double d_minus_i = ri - i / (double)n;

            d_plus = fmax(d_plus, d_plus_i);
            d_minus = fmax(d_minus, d_minus_i);

            printf("| %4d | %6.2f | %8.2f | %8.2f | %15.2f | %15.2f |\n", i + 1, ri, (i + 1.0) / n, i / (double)n, d_plus_i, d_minus_i);
        }

        printf("D+ = %f, D- = %f\n", d_plus, d_minus);
        return fmax(d_plus, d_minus);
    }

    // Function to get the critical value from the KS table
    double ks_critical_value(int n, double alpha)
    {
        if (alpha == 0.1)
            return 1.22 / sqrt(n);
        if (alpha == 0.05)
            return 1.36 / sqrt(n);
        if (alpha == 0.01)
            return 1.63 / sqrt(n);
        return -1.0; // Invalid alpha
    }

    int main()
    {
        int n;
        double alpha;

        // Get the sample size
        printf("Enter the sample size: ");
        scanf("%d", &n);

        // Allocate memory for the sample
        double *sample = malloc(n * sizeof(double));

        // Get the sample values
        printf("Enter the sample values:\n");
        for (int i = 0; i < n; i++)
        {
            printf("Value %d: ", i + 1);
            scanf("%lf", &sample[i]);
        }

        // Get the significance level
        printf("Enter the significance level (e.g., 0.05): ");
        scanf("%lf", &alpha);

        // Perform the KS test
        double d = ks_test(sample, n);
        double d_alpha = ks_critical_value(n, alpha);

        // Print the results
        printf("D calculated: %f\n", d);
        printf("D critical (alpha = %f): %f\n", alpha, d_alpha);

        // Decision rule
        if (d < d_alpha)
        {
            printf("Do not reject the null hypothesis (H0)\n");
        }
        else
        {
            printf("Reject the null hypothesis (H0)\n");
        }

        // Free the allocated memory
        free(sample);

        return 0;
    }

    */

/*
// Lab-12
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_GAP 10 // Define the maximum gap length

// Function to calculate the theoretical CDF
double theoretical_cdf(double x)
{
    return 1 - pow(0.9, x + 1);
}

// Function to calculate the empirical CDF
double empirical_cdf(int *gaps, int size, int x)
{
    int count = 0;
    for (int i = 0; i < size; i++)
    {
        if (gaps[i] <= x)
            count++;
    }
    return size > 0 ? (double)count / size : 0;
}

int main()
{
    int n, target, max_gap;
    double alpha, D_alpha;

    // Input sequence size
    printf("Enter the number of elements in the sequence: ");
    scanf("%d", &n);

    int sequence[n];
    int gaps[n];
    int gap_counts[MAX_GAP + 1] = {0};

    // Input sequence values
    printf("Enter the sequence values: ");
    for (int i = 0; i < n; i++)
    {
        scanf("%d", &sequence[i]);
    }

    // Input target value
    printf("Enter the target value: ");
    scanf("%d", &target);

    // Input maximum gap to count
    printf("Enter the maximum gap to count: ");
    scanf("%d", &max_gap);

    // Input significance level alpha
    printf("Enter the significance level alpha: ");
    scanf("%lf", &alpha);

    // Input critical value D_alpha
    printf("Enter the critical value D_alpha: ");
    scanf("%lf", &D_alpha);

    int gap = 0, gap_index = 0;

    // Perform gap test
    for (int i = 0; i < n; i++)
    {
        if (sequence[i] == target)
        {
            if (gap > 0)
            {
                if (gap <= max_gap)
                {
                    gap_counts[gap]++;
                }
                else
                {
                    gap_counts[max_gap]++;
                }
                gaps[gap_index++] = gap;
            }
            gap = 0;
        }
        else
        {
            gap++;
        }
    }

    // Calculate and print D
    double D = 0.0;
    printf("\nCalculations:\n x |   F(x)   |  S_N(x)  | |F(x) - S_N(x)|\n---|----------|----------|---------------\n");
    for (int i = 0; i <= max_gap; i++)
    {
        double F_x = theoretical_cdf(i);
        double S_N_x = empirical_cdf(gaps, gap_index, i);
        double D_current = fabs(F_x - S_N_x);
        printf("%2d   %8.6lf   %8.6lf   %13.6lf\n", i, F_x, S_N_x, D_current);
        if (D_current > D)
            D = D_current;
    }

    // Output result
    printf("\nMaximum Deviation D: %lf\n", D);
    printf("D_cal %s D_alpha, Null hypothesis is %s.\n", (D >= D_alpha) ? ">=" : "<", (D >= D_alpha) ? "rejected" : "not rejected");

    return 0;
}

*/