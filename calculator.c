#include <complex.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_INPUT_LEN 1024
#define PI 3.14159265358979323846L

typedef double complex cplx_t;

typedef struct {
  char *integer;
  char *fraction;
  bool is_negative;
} decimal_t;

decimal_t *decimal_construct(const char *integer, const char *fraction,
                             bool is_negative) {
  decimal_t *copy = malloc(sizeof(decimal_t));
  if (!copy) return NULL;

  copy->integer = (integer != NULL) ? strdup(integer) : "";
  if (!copy->integer) {
    free(copy);
    return NULL;
  }

  copy->fraction = (fraction != NULL) ? strdup(fraction) : "";
  if (!copy->fraction) {
    free(copy);
    return NULL;
  }

  copy->is_negative = is_negative;
  return copy;
}

void calculate(int argc, char *argv[]);
void interactive_mode();
void reverse(char *str);
char *convert_scientific(const char *input);
char *format_scientific(const char *num, bool approximate);
decimal_t *parse_number(const char *input);
void free_decimal(decimal_t *num);
char *add_numbers(decimal_t *a, decimal_t *b);
char *subtract_numbers(decimal_t *a, decimal_t *b);
char *multiply_numbers(decimal_t *a, decimal_t *b);
char *divide_numbers(decimal_t *a, decimal_t *b);
char *handle_sign(char *num, bool is_negative);
char *remove_leading_zeros(char *str);
char *remove_trailing_zeros(char *str);
void align_decimals(decimal_t *a, decimal_t *b);
char *add_strings(const char *a, const char *b, bool a_neg, bool b_neg,
                  bool *result_neg);
char *subtract_strings(const char *a, const char *b);
char *fft_multiply(const char *a, const char *b);
int compare_abs(const char *a, const char *b);
char *add_decimal_point(const char *num, int decimal_pos);
int count_leading_zeros(const char *num);
char *adjust_precision(const char *num, int precision);
char *reciprocal_newton(decimal_t *num);
decimal_t *change_format(char *input);
char *square_root(decimal_t *num);

// Data operation
decimal_t *parse_number(const char *input) {
  char *converted = convert_scientific(input);
  if (converted == NULL) {
    free(converted);
    return decimal_construct(NULL, NULL, false);
  }
  decimal_t *num = decimal_construct(NULL, NULL, false);
  char buffer[MAX_INPUT_LEN];
  strncpy(buffer, converted, MAX_INPUT_LEN - 1);
  buffer[MAX_INPUT_LEN - 1] = '\0';
  free(converted);

  if (buffer[0] == '-') {
    num->is_negative = true;
    memmove(buffer, buffer + 1, strlen(buffer) + 1);
  }

  char *p = buffer;
  int dot_count = 0;
  bool valid = true;

  if (*p == '\0') valid = false;

  char *start = p;
  while (*p && valid) {
    if (*p == '.') {
      if (++dot_count > 1 || p == start || *(p + 1) == '\0') {
        valid = false;
      }
    } else if (!isdigit(*p)) {
      valid = false;
    }
    ++p;
  }

  if (!valid) {
    num->integer = strdup("0");
    num->fraction = strdup("0");
    return num;
  }

  char *dot = strchr(buffer, '.');
  if (dot) {
    *dot = '\0';
    num->integer = strdup(buffer[0] ? buffer : "0");
    num->fraction = strdup(dot + 1);
  } else {
    num->integer = strdup(buffer);
    num->fraction = strdup("");
  }

  num->integer = remove_leading_zeros(num->integer);
  num->fraction = remove_trailing_zeros(num->fraction);
  return num;
}
char *convert_scientific(const char *input) {
  char *str = strdup(input);
  if (!str) return NULL;

  char *e_ptr = strpbrk(str, "eE");
  if (!e_ptr) return str;

  *e_ptr = '\0';
  char *mantissa = str;
  char *exp_str = e_ptr + 1;

  if (*exp_str == '\0' ||
      (!isdigit(*exp_str) && *exp_str != '+' && *exp_str != '-')) {
    free(str);
    return NULL;
  }

  bool exp_neg = false;
  if (*exp_str == '+')
    exp_str++;
  else if (*exp_str == '-') {
    exp_neg = true;
    exp_str++;
  }

  char *endptr;
  long exponent = strtol(exp_str, &endptr, 10);
  if (*endptr != '\0') {
    free(str);
    return NULL;
  }
  char *dot = strchr(mantissa, '.');
  int original_len = strlen(mantissa);
  int decimal_pos = dot ? (dot - mantissa + 1) : original_len;

  if (dot) {
    memmove(dot, dot + 1, strlen(dot));
    decimal_pos--;
  }

  int new_decimal_pos = decimal_pos + (exp_neg ? -exponent : exponent);
  int total_digits = strlen(mantissa);

  char *result = calloc(total_digits + abs(new_decimal_pos) + 3, 1);
  if (!result) {
    free(str);
    return NULL;
  }

  if (new_decimal_pos > 0) {
    if (total_digits < new_decimal_pos) {
      strcpy(result, mantissa);
      memset(result + total_digits, '0', new_decimal_pos - total_digits);
    } else {
      strncpy(result, mantissa, new_decimal_pos);
      if (total_digits > new_decimal_pos) {
        result[new_decimal_pos] = '.';
        strcpy(result + new_decimal_pos + 1, mantissa + new_decimal_pos);
      }
    }
    result = remove_leading_zeros(result);
  } else {
    strcpy(result, "0.");
    if (new_decimal_pos < 0) {
      int zeros_to_add = -new_decimal_pos;
      if (zeros_to_add > 0) {
        memset(result + 2, '0', zeros_to_add);
      }
      strcat(result + 2 + zeros_to_add, mantissa);
    } else {
      strcat(result, mantissa);
    }
  }
  free(str);
  return result;
}

char *format_scientific(const char *num, bool approximate) {
  if (!num || !*num) return strdup("0");

  char *processed = strdup(num);
  if (!processed) return NULL;

  bool is_neg = (processed[0] == '-');
  if (is_neg) memmove(processed, processed + 1, strlen(processed));
  processed = remove_leading_zeros(processed);
  if (processed[0] == '.') {
    char *temp = malloc(strlen(processed) + 20);
    if (!temp) {
      free(processed);
      return NULL;
    }
    strcpy(temp, "0");
    strcat(temp, processed);
    free(processed);
    processed = temp;
  } else if (strcmp(processed, "") == 0) {
    free(processed);
    processed = strdup("0");
  }
  char *first_digit = processed;
  while (*first_digit == '0' || *first_digit == '.') {
    if (*first_digit == '\0') break;
    first_digit++;
  }

  if (*first_digit == '\0') {
    free(processed);
    return strdup("0");
  }

  char *dot = strchr(processed, '.');
  int exp = -1;
  bool use_sci = false;
  bool has_carry = false;
  bool shorten = false;
  if (approximate && strlen(processed) > 10) {
    shorten = true;
  }
  if (dot != NULL) {
    if (processed[0] == '0' || processed[0] == '.') {
      exp = strlen(first_digit) - strlen(dot);
    } else {
      exp = strlen(first_digit) - strlen(dot) - 1;
    }
  } else {
    exp = strlen(first_digit) - 1;
  }
  use_sci = (exp >= 10 || exp <= -10);
  if (!use_sci) {
    if (shorten) {
      int new_len = strlen(processed) - strlen(dot) + 10;
      char *temp = malloc(new_len + 1);
      strncpy(temp, processed, new_len);
      temp[new_len] = '\0';
      if (processed[new_len] >= '5') {
        has_carry = true;
      }
      char *processed_temp = processed;
      processed = temp;
      free(processed_temp);
      dot = strchr(processed, '.');
      first_digit = processed;
    }
    char *res = malloc(strlen(processed) + is_neg + 1);
    if (!res) {
      free(processed);
      return NULL;
    }
    sprintf(res, "%s%s", is_neg ? "-" : "", processed);
    free(processed);
    if (has_carry) {
      decimal_t *num_carry = decimal_construct("", "000000001", false);
      decimal_t *res_num = parse_number(res);
      char *temp = add_numbers(num_carry, res_num);
      free(res);
      free_decimal(res_num);
      free_decimal(num_carry);
      res = temp;
    }
    return res;
  }
  if (*first_digit == '.') {
    first_digit++;
  }
  first_digit = remove_leading_zeros(first_digit);
  char *mantissa = calloc(strlen(processed) + 2, 1);
  int pos = 0;
  mantissa[pos++] = *first_digit;
  if (*(first_digit + 1) != '\0') {
    mantissa[pos++] = '.';
    char *p = first_digit + 1;
    while (*p) {
      if (*p != '.') mantissa[pos++] = *p;
      p++;
    }
  }
  mantissa[pos] = '\0';
  remove_trailing_zeros(mantissa);
  if (mantissa[pos - 1] == '.') mantissa[pos - 1] = '\0';

  char *result = malloc(strlen(mantissa) + 200);
  if (shorten) {
    int new_len = 11;
    char *temp = malloc(new_len + 1);
    strncpy(temp, mantissa, new_len);
    temp[new_len] = '\0';
    if (mantissa[new_len] >= '5') {
      has_carry = true;
    }
    char *mantissa_temp = mantissa;
    mantissa = temp;
    free(mantissa_temp);
  }
  if (has_carry) {
    decimal_t *num_carry = decimal_construct("", "000000001", false);
    decimal_t *res_num = parse_number(mantissa);
    char *temp = add_numbers(num_carry, res_num);
    if (temp[2] == '.') {
      temp[1] = '.';
      temp[2] = '0';
      exp++;
    }
    free(mantissa);
    free_decimal(res_num);
    free_decimal(num_carry);
    mantissa = temp;
  }
  sprintf(result, "%s%se%02d", is_neg ? "-" : "", mantissa, exp);
  free(processed);
  free(mantissa);
  return result;
}

// Main function
void calculate(int argc, char *argv[]) {
  if (strcmp(argv[1], "sqrt") == 0) {
    if (argv[2] == NULL) {
      fprintf(stderr, "Error: Invalid expression format\n");
      return;
    }
    if (argv[3] != NULL) {
      fprintf(stderr, "Error: Invalid expression format\n");
      return;
    }
    if (argv[2][0] == '-') {
      fprintf(stderr, "Error: The input should NOT be negative!\n");
      return;
    }
    bool exit_tag = false;
    for (int i = 0; i < strlen(argv[2]); i++) {
      if ((argv[2][i] < '0' || argv[2][i] > '9') && argv[2][i] != '.' &&
          argv[2][i] != 'e' && argv[2][i] != 'E' && argv[2][i] != '-') {
        exit_tag = true;
        continue;
      }
    }
    if (exit_tag) {
      fprintf(stderr, "Error: Invalid expression format\n");
      return;
    }
    decimal_t *num = parse_number(argv[2]);
    int len_num = strlen(num->integer) + strlen(num->fraction);
    int len_num_fraction = strlen(num->fraction);
    if (len_num > 500) {
      fprintf(stderr,
              "Error: The input number has too many significant digits. Please "
              "enter a valid number.\n");
      free_decimal(num);
      return;
    }
    if (len_num_fraction > 18) {
      fprintf(stderr,
              "Error: The input number has too many significant digits. Please "
              "enter a valid number.\n");
      free_decimal(num);
      return;
    }
    char *result = square_root(num);
    if (result) {
      char *sci_result = format_scientific(result, true);
      printf("Result: %s\n", sci_result);
      free(sci_result);
      free(result);
    }
    free_decimal(num);
    return;
  }
  for (int i = 0; i < strlen(argv[1]); i++) {
    if ((argv[1][i] < '0' || argv[1][i] > '9') && argv[1][i] != '.' &&
        argv[1][i] != 'e' && argv[1][i] != 'E' && argv[1][i] != '-') {
      fprintf(stderr,
              "Error: Invalid expression format. Usage: %s [num1] [operator] "
              "[num2]\n",
              argv[0]);
      return;
    }
  }
  for (int i = 0; i < strlen(argv[3]); i++) {
    if ((argv[3][i] < '0' || argv[3][i] > '9') && argv[3][i] != '.' &&
        argv[3][i] != 'e' && argv[3][i] != 'E' && argv[3][i] != '-') {
      fprintf(stderr,
              "Error: Invalid expression format. Usage: %s [num1] [operator] "
              "[num2]\n",
              argv[0]);
      return;
    }
  }
  decimal_t *num1 = parse_number(argv[1]);
  decimal_t *num2 = parse_number(argv[3]);
  const char *op = argv[2];

  int len1 = strlen(num1->integer) + strlen(num1->fraction);
  int len2 = strlen(num2->integer) + strlen(num2->fraction);
  int len2_integer = strlen(num2->integer);
  int len2_fraction = strlen(num2->fraction);
  if (len1 > 1000 || len2 > 1000) {
    fprintf(stderr,
            "Error: The input number has too many significant digits. Please "
            "enter a valid number.\n");
    free_decimal(num1);
    free_decimal(num2);
    return;
  }
  if (argv[2][0] == '/' && len2_integer > 500) {
    fprintf(stderr,
            "Error: The input number has too many significant digits. Please "
            "enter a valid number.\n");
    free_decimal(num1);
    free_decimal(num2);
    return;
  }
  if (argv[2][0] == '/' && len2_fraction > 500) {
    fprintf(stderr,
            "Error: The input number has too many significant digits. Please "
            "enter a valid number.\n");
    free_decimal(num1);
    free_decimal(num2);
    return;
  }
  if (num1->integer == NULL || num2->integer == NULL) {
    fprintf(stderr, "Error: Invalid scientific notation format\n");
    free_decimal(num1);
    free_decimal(num2);
    return;
  }

  char *result = NULL;
  if (strcmp(op, "+") == 0) {
    result = add_numbers(num1, num2);
  } else if (strcmp(op, "-") == 0) {
    result = subtract_numbers(num1, num2);
  } else if (strcmp(op, "x") == 0) {
    result = multiply_numbers(num1, num2);
  } else if (strcmp(op, "/") == 0) {
    result = divide_numbers(num1, num2);
  } else {
    fprintf(stderr, "Error: Invalid operator\n");
  }

  if (result) {
    char *sci_result = format_scientific(result, true);
    printf("%s %s %s = %s\n", argv[1], op, argv[3], sci_result);
    free(sci_result);
    free(result);
  }

  free_decimal(num1);
  free_decimal(num2);
  return;
}
void interactive_mode() {
  printf("High Precision Calculator (enter 'quit' to exit)\n");
  printf(
      "Supports +, -, x, / operations and scientific notation (e.g. 1.2e3)\n");
  printf(
      "Supports output mode precise and approximate, use command [mode "
      "precise] and [mode approximate] to shift\n");
  char input[MAX_INPUT_LEN];
  bool approximate = true;
  while (1) {
    printf("> ");
    if (!fgets(input, MAX_INPUT_LEN, stdin)) break;
    input[strcspn(input, "\n")] = '\0';

    if (strcmp(input, "quit") == 0) break;

    char *tokens[3] = {NULL, NULL, NULL};
    char *token = strtok(input, " ");
    int i = 0;
    while (token != NULL && i < 3) {
      tokens[i++] = token;
      token = strtok(NULL, " ");
    }
    if (token != NULL) {
      fprintf(stderr, "Error: Invalid expression format\n");
      continue;
    }
    if (tokens[0] == NULL) {
      continue;
    }
    if (strcmp(tokens[0], "sqrt") == 0) {
      if (tokens[1] == NULL) {
        fprintf(stderr, "Error: Invalid expression format\n");
        continue;
      }
      if (tokens[2] != NULL) {
        fprintf(stderr, "Error: Invalid expression format\n");
        continue;
      }
      bool exit_tag = false;
      if (tokens[1][0] == '-') {
        fprintf(stderr, "Error: The input should NOT be negative!\n");
        continue;
      }
      for (int i = 0; i < strlen(tokens[1]); i++) {
        if ((tokens[1][i] < '0' || tokens[1][i] > '9') && tokens[1][i] != '.' &&
            tokens[1][i] != 'e' && tokens[1][i] != 'E' && tokens[1][i] != '-') {
          exit_tag = true;
          continue;
        }
      }
      if (exit_tag) {
        fprintf(stderr, "Error: Invalid expression format\n");
        continue;
      }
      decimal_t *num = parse_number(tokens[1]);
      int len_num = strlen(num->integer) + strlen(num->fraction);
      int len_num_fraction = strlen(num->fraction);
      if (len_num > 500) {
        fprintf(stderr,
                "Error: The input number has too many significant digits. "
                "Please enter a valid number.\n");
        free_decimal(num);
        continue;
      }
      if (len_num_fraction > 18) {
        fprintf(stderr,
                "Error: The input number has too many significant digits. "
                "Please enter a valid number.\n");
        free_decimal(num);
        continue;
      }
      char *result = square_root(num);
      if (result) {
        char *sci_result = format_scientific(result, approximate);
        printf("Result: %s\n", sci_result);
        free(sci_result);
        free(result);
      }
      free_decimal(num);
      continue;
    }
    if (strcmp(tokens[0], "mode") == 0) {
      if (tokens[1] == NULL) {
        printf(
            "Error: Mode swift expression format: [mode precise] or [mode "
            "approximate]\n");
      } else if ((strcmp(tokens[1], "precise") != 0 &&
                  strcmp(tokens[1], "approximate") != 0) ||
                 tokens[2] != NULL) {
        printf(
            "Error: Mode swift expression format: [mode precise] or [mode "
            "approximate]\n");
      } else if (strcmp(tokens[1], "precise") == 0) {
        approximate = false;
        printf("Switch to precise mode\n");
      } else if (strcmp(tokens[1], "approximate") == 0) {
        approximate = true;
        printf("Switch to approximate mode\n");
      }
      continue;
    }
    if (tokens[2] == NULL) {
      fprintf(stderr, "Error: Invalid expression format\n");
      continue;
    }
    bool invalid_input = false;
    for (int i = 0; i < strlen(tokens[0]); i++) {
      if ((tokens[0][i] < '0' || tokens[0][i] > '9') && tokens[0][i] != '.' &&
          tokens[0][i] != 'e' && tokens[0][i] != 'E' && tokens[0][i] != '-') {
        invalid_input = true;
        continue;
      }
    }
    for (int i = 0; i < strlen(tokens[2]); i++) {
      if ((tokens[2][i] < '0' || tokens[2][i] > '9') && tokens[2][i] != '.' &&
          tokens[2][i] != 'e' && tokens[2][i] != 'E' && tokens[2][i] != '-') {
        invalid_input = true;
        continue;
      }
    }
    if (invalid_input) {
      fprintf(stderr, "Error: Invalid expression format\n");
      continue;
    }

    if (i != 3 || strlen(tokens[1]) != 1 || !strchr("+-x/", tokens[1][0])) {
      fprintf(stderr,
              "Error: Invalid expression format. Format: [num1] [operator] "
              "[num2]\n");
      continue;
    }

    decimal_t *num1 = parse_number(tokens[0]);
    decimal_t *num2 = parse_number(tokens[2]);

    int len1 = strlen(num1->integer) + strlen(num1->fraction);
    int len2 = strlen(num2->integer) + strlen(num2->fraction);
    int len2_integer = strlen(num2->integer);
    int len2_fraction = strlen(num2->fraction);
    if (len1 > 1000 || len2 > 1000) {
      fprintf(stderr,
              "Error: The input number has too many significant digits. Please "
              "enter a valid number.\n");
      free_decimal(num1);
      free_decimal(num2);
      continue;
    }
    if (tokens[1][0] == '/' && len2_integer > 500) {
      fprintf(stderr,
              "Error: The input number has too many significant digits. Please "
              "enter a valid number.\n");
      free_decimal(num1);
      free_decimal(num2);
      continue;
    }
    if (tokens[1][0] == '/' && len2_fraction > 500) {
      fprintf(stderr,
              "Error: The input number has too many significant digits. Please "
              "enter a valid number.\n");
      free_decimal(num1);
      free_decimal(num2);
      continue;
    }
    if (num1->integer == NULL || num2->integer == NULL) {
      fprintf(stderr, "Error: Invalid scientific notation format\n");
      free_decimal(num1);
      free_decimal(num2);
      continue;
    }

    char *result = NULL;
    switch (tokens[1][0]) {
      case '+':
        result = add_numbers(num1, num2);
        break;
      case '-':
        result = subtract_numbers(num1, num2);
        break;
      case 'x':
      case '*':
        result = multiply_numbers(num1, num2);
        break;
      case '/':
        result = divide_numbers(num1, num2);
        break;
      default:
        fprintf(stderr, "Error: Unsupported operator\n");
    }

    if (result) {
      char *sci_result = format_scientific(result, approximate);
      printf("Result: %s\n", sci_result);
      free(sci_result);
      free(result);
    }

    free_decimal(num1);
    free_decimal(num2);
  }
  return;
}
int main(int argc, char *argv[]) {
  if (argc == 4) {
    calculate(argc, argv);
  } else if (argc == 3 && strcmp(argv[1], "sqrt") == 0) {
    calculate(argc, argv);
  } else if (argc != 1) {
    fprintf(stderr, "Usage: %s [num1] [operator] [num2]\n", argv[0]);
    return 1;
  } else {
    interactive_mode();
  }
  return 0;
}

// Calculate operation
char *add_numbers(decimal_t *a, decimal_t *b) {
  align_decimals(a, b);
  const int frac_len = strlen(a->fraction);

  char *num1 = malloc(strlen(a->integer) + frac_len + 1);
  sprintf(num1, "%s%s", a->integer, a->fraction);

  char *num2 = malloc(strlen(b->integer) + frac_len + 1);
  sprintf(num2, "%s%s", b->integer, b->fraction);

  bool result_neg;
  char *abs_res =
      add_strings(num1, num2, a->is_negative, b->is_negative, &result_neg);
  free(num1);
  free(num2);

  int total_len = strlen(abs_res);
  int int_len = total_len - frac_len;

  if (int_len < 0) {
    int needed_zeros = -int_len;
    char *padded = malloc(frac_len + 1);
    memset(padded, '0', needed_zeros);
    strcpy(padded + needed_zeros, abs_res);
    free(abs_res);
    abs_res = padded;
    total_len = strlen(abs_res);
    int_len = total_len - frac_len;
  }

  char *result = malloc(total_len + 2);
  strncpy(result, abs_res, int_len);
  result[int_len] = '\0';

  if (frac_len > 0) {
    result[int_len] = '.';
    strncpy(result + int_len + 1, abs_res + int_len, frac_len);
    result[int_len + 1 + frac_len] = '\0';
  }

  free(abs_res);
  result = remove_leading_zeros(result);

  return handle_sign(result, result_neg);
}

char *subtract_numbers(decimal_t *a, decimal_t *b) {
  decimal_t *temp = decimal_construct(b->integer, b->fraction, !b->is_negative);
  char *result = add_numbers(a, temp);
  free_decimal(temp);
  return result;
}

char *multiply_numbers(decimal_t *a, decimal_t *b) {
  char *num1 = malloc(strlen(a->integer) + strlen(a->fraction) + 1);
  sprintf(num1, "%s%s", a->integer, a->fraction);
  char *num2 = malloc(strlen(b->integer) + strlen(b->fraction) + 1);
  sprintf(num2, "%s%s", b->integer, b->fraction);

  const int decimal_places = strlen(a->fraction) + strlen(b->fraction);
  char *product = fft_multiply(num1, num2);
  free(num1);
  free(num2);

  int total_len = strlen(product);
  int int_len = total_len - decimal_places;

  if (int_len < 0) {
    int leading_zeros = -int_len;
    char *temp = malloc(total_len + leading_zeros + 1);
    memset(temp, '0', leading_zeros);
    strcpy(temp + leading_zeros, product);
    free(product);
    product = temp;
    total_len = strlen(product);
    int_len = 0;
  }

  char *result = malloc(total_len + 2);
  strncpy(result, product, int_len);
  result[int_len] = '\0';

  if (decimal_places > 0) {
    if (int_len == 0) {
      result = realloc(result, decimal_places + 3);
      strcpy(result, "0.");
      strncpy(result + 2, product, decimal_places);
      result[2 + decimal_places] = '\0';
    } else {
      result[int_len] = '.';
      strncpy(result + int_len + 1, product + int_len, decimal_places);
      result[int_len + 1 + decimal_places] = '\0';
    }
    remove_trailing_zeros(result + (int_len == 0 ? 2 : int_len + 1));

    if (result[int_len + 1] == '\0') {
      result[int_len] = '\0';
    }
  }

  free(product);
  return handle_sign(remove_leading_zeros(result),
                     a->is_negative != b->is_negative);
}

void fft(cplx_t a[], int n, bool invert) {
  if (n == 1) return;

  cplx_t *a0 = malloc(n / 2 * sizeof(cplx_t));
  cplx_t *a1 = malloc(n / 2 * sizeof(cplx_t));

  for (int i = 0; i < n / 2; ++i) {
    a0[i] = a[2 * i];
    a1[i] = a[2 * i + 1];
  }

  fft(a0, n / 2, invert);
  fft(a1, n / 2, invert);

  const double angle = 2 * PI / n * (invert ? -1 : 1);
  cplx_t w = 1.0;
  const cplx_t wn = cexp(angle * I);

  for (int i = 0; i < n / 2; ++i) {
    const cplx_t t = w * a1[i];
    a[i] = a0[i] + t;
    a[i + n / 2] = a0[i] - t;
    if (invert) {
      a[i] /= 2;
      a[i + n / 2] /= 2;
    }
    w *= wn;
  }

  free(a0);
  free(a1);
}

char *fft_multiply(const char *a, const char *b) {
  const int len1 = strlen(a);
  const int len2 = strlen(b);
  int n = 1;
  while (n < len1 + len2) n <<= 1;

  cplx_t *fa = calloc(n, sizeof(cplx_t));
  cplx_t *fb = calloc(n, sizeof(cplx_t));

  for (int i = 0; i < len1; ++i) fa[i] = a[len1 - 1 - i] - '0';
  for (int i = 0; i < len2; ++i) fb[i] = b[len2 - 1 - i] - '0';

  fft(fa, n, false);
  fft(fb, n, false);

  for (int i = 0; i < n; ++i) fa[i] *= fb[i];

  fft(fa, n, true);

  int *res = calloc(n, sizeof(int));
  for (int i = 0; i < n; ++i) {
    res[i] = (int)(creal(fa[i]) + 0.5);
  }

  for (int i = 0; i < n - 1; ++i) {
    res[i + 1] += res[i] / 10;
    res[i] %= 10;
  }

  int pos = n - 1;
  while (pos > 0 && res[pos] == 0) --pos;

  char *result = malloc(pos + 2);
  for (int i = 0; i <= pos; ++i) {
    result[i] = res[pos - i] + '0';
  }
  result[pos + 1] = '\0';

  free(fa);
  free(fb);
  free(res);
  return remove_leading_zeros(result);
}

char *divide_numbers(decimal_t *a, decimal_t *b) {
  if (strcmp(b->integer, "0") == 0 && strcmp(b->fraction, "") == 0) {
    fprintf(stderr, "Error: Division by zero\n");
    return NULL;
  }
  bool b_neg = b->is_negative;
  if (b_neg) {
    b->is_negative = false;
  }
  decimal_t *reciprocal = change_format(reciprocal_newton(b));
  char *product = multiply_numbers(a, reciprocal);
  free_decimal(reciprocal);
  return handle_sign(product, a->is_negative != b_neg);
}

char *reciprocal_newton(decimal_t *num) {
  decimal_t *x =
      decimal_construct(num->integer, num->fraction, num->is_negative);
  decimal_t *num_0_5 = decimal_construct("0", "5", false);
  decimal_t *num_0_96 = decimal_construct("0", "96", false);
  decimal_t *num_1 = decimal_construct("1", "", false);
  decimal_t *num_1_48 = decimal_construct("1", "48", false);
  decimal_t *num_2 = decimal_construct("2", "", false);
  int scale = 0;
  bool need_scale = true;
  while (need_scale) {
    decimal_t *cmp_half = change_format(subtract_numbers(x, num_0_5));
    decimal_t *cmp_one = change_format(subtract_numbers(x, num_1));
    if (cmp_half->is_negative) {
      decimal_t *x_current = x;
      x = change_format(multiply_numbers(x, num_2));
      scale++;
      free_decimal(x_current);
    } else if (!cmp_one->is_negative) {
      decimal_t *x_current = x;
      x = change_format(multiply_numbers(x, num_0_5));
      scale--;
      free_decimal(x_current);
    } else {
      need_scale = false;
    }
    free_decimal(cmp_half);
    free_decimal(cmp_one);
  }
  decimal_t *zero_96_x = change_format(multiply_numbers(x, num_0_96));
  decimal_t *x_current = x;
  x = change_format(subtract_numbers(num_1_48, zero_96_x));
  free_decimal(x_current);
  free_decimal(zero_96_x);

  if (scale < 0) {
    for (int i = 0; i < -scale; i++) {
      decimal_t *x_current = x;
      x = change_format(multiply_numbers(x, num_0_5));
      free_decimal(x_current);
    }
  } else if (scale > 0) {
    for (int i = 0; i < scale; i++) {
      decimal_t *x_current = x;
      x = change_format(multiply_numbers(x, num_2));
      free_decimal(x_current);
    }
  }

  for (int i = 0; i < 14; i++) {
    decimal_t *x_squared = change_format(multiply_numbers(x, x));
    decimal_t *x_squared_n = change_format(multiply_numbers(x_squared, num));
    decimal_t *two_x = change_format(multiply_numbers(x, num_2));
    decimal_t *x_current = x;
    x = change_format(subtract_numbers(two_x, x_squared_n));
    free_decimal(x_current);

    free_decimal(x_squared);
    free_decimal(x_squared_n);
    free_decimal(two_x);
  }
  char *result = multiply_numbers(x, num_1);

  free_decimal(num_0_5);
  free_decimal(num_0_96);
  free_decimal(num_1);
  free_decimal(num_1_48);
  free_decimal(num_2);
  free_decimal(x);
  return result;
}
char *square_root(decimal_t *num) {
  decimal_t *simulated_num =
      decimal_construct(num->integer, num->fraction, num->is_negative);
  decimal_t *num_0_5 = decimal_construct("0", "5", false);
  decimal_t *num_0_7 = decimal_construct("0", "7", false);
  decimal_t *num_1 = decimal_construct("1", "0", false);
  decimal_t *num_1_5 = decimal_construct("1", "5", false);
  decimal_t *num_2 = decimal_construct("2", "0", false);
  decimal_t *square_num =
      change_format(multiply_numbers(simulated_num, simulated_num));

  decimal_t *cmp_square = change_format(subtract_numbers(num, square_num));
  if (!cmp_square->is_negative) {
    while (!cmp_square->is_negative) {
      decimal_t *simulated_num_now = simulated_num;
      simulated_num =
          change_format(multiply_numbers(simulated_num_now, num_1_5));
      free_decimal(simulated_num_now);
      decimal_t *square_num_now = square_num;
      square_num =
          change_format(multiply_numbers(simulated_num, simulated_num));
      free_decimal(square_num_now);
      decimal_t *cmp_square_now = cmp_square;
      cmp_square = change_format(subtract_numbers(num, square_num));
      free_decimal(cmp_square_now);
    }
  } else {
    while (cmp_square->is_negative) {
      decimal_t *simulated_num_now = simulated_num;
      simulated_num =
          change_format(multiply_numbers(simulated_num_now, num_0_7));
      free_decimal(simulated_num_now);
      decimal_t *square_num_now = square_num;
      square_num =
          change_format(multiply_numbers(simulated_num, simulated_num));
      free_decimal(square_num_now);
      decimal_t *cmp_square_now = cmp_square;
      cmp_square = change_format(subtract_numbers(num, square_num));
      free_decimal(cmp_square_now);
    }
  }
  decimal_t *num_reciprocal =
      change_format(divide_numbers(num_1, simulated_num));
  free_decimal(cmp_square);
  free_decimal(square_num);

  for (int i = 0; i < 14; i++) {
    decimal_t *square_2_num =
        change_format(multiply_numbers(num_reciprocal, num_reciprocal));
    decimal_t *square_3_num =
        change_format(multiply_numbers(square_2_num, num_reciprocal));
    decimal_t *valued_num =
        change_format(multiply_numbers(num_reciprocal, num_1_5));
    decimal_t *valued_square =
        change_format(multiply_numbers(square_3_num, num_0_5));
    decimal_t *valued_square_num =
        change_format(multiply_numbers(valued_square, num));
    decimal_t *num_reciprocal_now = num_reciprocal;
    char *temp = subtract_numbers(valued_num, valued_square_num);
    num_reciprocal = change_format(temp);
    free_decimal(num_reciprocal_now);
    free_decimal(square_2_num);
    free_decimal(square_3_num);
    free_decimal(valued_num);
    free_decimal(valued_square);
    free_decimal(valued_square_num);
  }
  char *result = multiply_numbers(num, num_reciprocal);
  free_decimal(num_reciprocal);
  free_decimal(simulated_num);
  free_decimal(num_0_5);
  free_decimal(num_0_7);
  free_decimal(num_1);
  free_decimal(num_1_5);
  free_decimal(num_2);
  return handle_sign(result, false);
}

// Tool function

char *add_decimal_point(const char *num, int decimal_pos) {
  int len = strlen(num);
  char *res = NULL;

  if (decimal_pos <= 0) {
    int zeros_needed = -decimal_pos;
    res = malloc(len + zeros_needed + 3);
    sprintf(res, "0.");
    for (int i = 0; i < zeros_needed; i++) strcat(res, "0");
    strcat(res, num);
  } else if (decimal_pos >= len) {
    res = malloc(decimal_pos + 32);
    strcpy(res, num);
    for (int i = len; i < decimal_pos; i++) strcat(res, "0");
  } else {
    res = malloc(len + 2);
    strncpy(res, num, len - decimal_pos);
    res[len - decimal_pos] = '.';
    strncpy(res + len - decimal_pos + 1, num + len - decimal_pos, decimal_pos);
    res[len + 1] = '\0';
  }
  return remove_trailing_zeros(res);
}

int count_leading_zeros(const char *num) {
  int count = 0;
  while (num[count] == '0') count++;
  return count;
}

char *adjust_precision(const char *num, int precision) {
  int current_len = strlen(num);
  char *res = malloc(precision + 1);
  if (!res) return NULL;

  if (current_len >= precision) {
    strncpy(res, num, precision);
  } else {
    strcpy(res, num);
    memset(res + current_len, '0', precision - current_len);
  }
  res[precision] = '\0';
  return res;
}

char *handle_sign(char *num, bool neg) {
  if (neg && num[0] != '-') {
    if (num[0] == '.') {
      char *res = malloc(strlen(num) + 30);
      sprintf(res, "-0%s", num);
      free(num);
      return res;
    }
    char *res = malloc(strlen(num) + 2);
    sprintf(res, "-%s", num);
    free(num);
    return res;
  }
  if (num[0] == '.') {
    char *res = malloc(strlen(num) + 30);
    sprintf(res, "0%s", num);
    free(num);
    return res;
  }
  return num;
}

char *remove_leading_zeros(char *str) {
  if (!str || *str == '\0') return str;
  char *p = str;
  while (*p == '0' && *(p + 1) != '\0' && *(p + 1) != '.') ++p;
  memmove(str, p, strlen(p) + 1);
  return str;
}

char *remove_trailing_zeros(char *str) {
  if (!str || *str == '\0') return str;
  char *end = str + strlen(str) - 1;
  while (end >= str && *end == '0') *end-- = '\0';
  if (end >= str && *end == '.') *end = '\0';
  return str;
}

void free_decimal(decimal_t *num) {
  if (num->integer) free(num->integer);
  if (num->fraction) free(num->fraction);
  free(num);
}

void reverse(char *str) {
  int len = strlen(str);
  for (int i = 0; i < len / 2; ++i) {
    char t = str[i];
    str[i] = str[len - i - 1];
    str[len - i - 1] = t;
  }
}

int compare_abs(const char *a, const char *b) {
  char *pa = strdup(a);
  char *pb = strdup(b);
  pa = remove_leading_zeros(pa);
  pb = remove_leading_zeros(pb);

  char *dot_a = strchr(pa, '.');
  char *dot_b = strchr(pb, '.');
  int int_a_len = dot_a ? (dot_a - pa) : strlen(pa);
  int int_b_len = dot_b ? (dot_b - pb) : strlen(pb);

  if (int_a_len != int_b_len) {
    free(pa);
    free(pb);
    return int_a_len - int_b_len;
  }

  for (int i = 0; i < int_a_len; i++) {
    if (pa[i] != pb[i]) {
      int res = pa[i] - pb[i];
      free(pa);
      free(pb);
      return res;
    }
  }

  char *frac_a = dot_a ? (dot_a + 1) : "";
  char *frac_b = dot_b ? (dot_b + 1) : "";
  int max_frac =
      strlen(frac_a) > strlen(frac_b) ? strlen(frac_a) : strlen(frac_b);

  for (int i = 0; i < max_frac; i++) {
    char ca = (i < strlen(frac_a)) ? frac_a[i] : '0';
    char cb = (i < strlen(frac_b)) ? frac_b[i] : '0';
    if (ca != cb) {
      free(pa);
      free(pb);
      return ca - cb;
    }
  }

  free(pa);
  free(pb);
  return 0;
}

void align_decimals(decimal_t *a, decimal_t *b) {
  const int a_len = strlen(a->fraction);
  const int b_len = strlen(b->fraction);
  const int max_len = (a_len > b_len) ? a_len : b_len;

  a->fraction = realloc(a->fraction, max_len + 1);
  b->fraction = realloc(b->fraction, max_len + 1);

  if (a_len < max_len) {
    memset(a->fraction + a_len, '0', max_len - a_len);
    a->fraction[max_len] = '\0';
  }
  if (b_len < max_len) {
    memset(b->fraction + b_len, '0', max_len - b_len);
    b->fraction[max_len] = '\0';
  }
}

char *add_strings(const char *a, const char *b, bool a_neg, bool b_neg,
                  bool *result_neg) {
  if (a_neg == b_neg) {
    *result_neg = a_neg;
    char *ra = strdup(a);
    reverse(ra);
    char *rb = strdup(b);
    reverse(rb);

    const int max_len = (strlen(a) > strlen(b)) ? strlen(a) : strlen(b);
    char *result = malloc(max_len + 2);
    int carry = 0, idx = 0;

    for (int i = 0; i < max_len || carry; ++i) {
      const int n1 = (i < strlen(ra)) ? (ra[i] - '0') : 0;
      const int n2 = (i < strlen(rb)) ? (rb[i] - '0') : 0;
      const int sum = n1 + n2 + carry;
      result[idx++] = (sum % 10) + '0';
      carry = sum / 10;
    }
    result[idx] = '\0';

    reverse(result);
    free(ra);
    free(rb);
    return remove_leading_zeros(result);
  }

  const int cmp = compare_abs(a, b);
  if (cmp == 0) return strdup("0");

  const char *larger = (cmp > 0) ? a : b;
  const char *smaller = (cmp > 0) ? b : a;
  *result_neg = (cmp > 0) ? a_neg : b_neg;

  return subtract_strings(larger, smaller);
}

char *subtract_strings(const char *a, const char *b) {
  char *ra = strdup(a);
  reverse(ra);
  char *rb = strdup(b);
  reverse(rb);

  const int len_a = strlen(ra);
  char *result = malloc(len_a + 2);
  int borrow = 0, idx = 0;

  for (int i = 0; i < len_a; ++i) {
    int digit_a = (ra[i] - '0') - borrow;
    borrow = 0;
    const int digit_b = (i < strlen(rb)) ? (rb[i] - '0') : 0;

    if (digit_a < digit_b) {
      digit_a += 10;
      borrow = 1;
    }
    result[idx++] = (digit_a - digit_b) + '0';
  }

  result[idx] = '\0';
  reverse(result);
  free(ra);
  free(rb);
  return remove_leading_zeros(result);
}
decimal_t *change_format(char *input) {
  char *adjusted = format_scientific(input, false);
  free(input);
  decimal_t *result = parse_number(adjusted);
  free(adjusted);
  return result;
}
