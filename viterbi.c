/*
 * Viterbi encode/decode in C.
 * https://zlearning.netlify.app/communication/ecc/eccviterbi
 *
 * qianfan Zhao <qianfanguijin@163.com>
 */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <getopt.h>

static int option_verbose = 0;

#define verbose_printf(level, fmt, ...)		do {	\
	if (option_verbose >= level)			\
		printf(fmt, ##__VA_ARGS__);		\
} while (0)

#define bitof(byte, n)			(((byte) >> (n)) & 1)

static size_t str01_to_buffer(const char *s, const char **endp, uint8_t *buf,
			      size_t bufsz, int lsb_first)
{
	size_t bytes = 0, binary_counts = 0;

	if (buf && bufsz)
		buf[0] = 0;

	for (; *s != '\0' && bytes < bufsz; s++) {
		uint8_t bit = binary_counts % 8;
		int n = *s - '0';

		if (n != 0 && n != 1)
			break;

		if (lsb_first)
			buf[bytes] |= (n << bit);
		else
			buf[bytes] |= (n << (7 - bit));

		binary_counts++;
		if ((binary_counts % 8) == 0) {
			bytes++;
			if (bytes < bufsz)
				buf[bytes] = 0;
		}
	}

	if (endp)
		*endp = s;

	return binary_counts;
}

struct fec_table {
	uint8_t		m;
	uint8_t		o;
	uint8_t		next_m;
};

#define FEC_TABLE(m, o, next_m)		{ m, o, next_m }
#define CONV_MAX_M			4

static const struct fec_table state_tables[2][CONV_MAX_M] = {
	{
		FEC_TABLE(0b00, 0b00, 0b00),
		FEC_TABLE(0b01, 0b11, 0b00),
		FEC_TABLE(0b10, 0b10, 0b01),
		FEC_TABLE(0b11, 0b01, 0b01),
	}, {
		FEC_TABLE(0b00, 0b11, 0b10),
		FEC_TABLE(0b01, 0b00, 0b10),
		FEC_TABLE(0b10, 0b01, 0b11),
		FEC_TABLE(0b11, 0b10, 0b11),
	}
};

static void foreach_bit(const uint8_t *buf, size_t buf_bits,
			void (*f)(int, void *), void *private_data)
{
	size_t byte_idx = 0, total_bits = 0;
	uint8_t b;

	for (size_t i = 0; i < buf_bits; i++) {
		uint8_t bit_idx = total_bits & 0b111;
		int bi;

		if (bit_idx == 0)
			b = buf[byte_idx++];

		bi = !!(b & (1 << (7 - bit_idx)));
		f(bi, private_data);
		total_bits++;
	}
}

static void print_bit(int bi, void *private_data)
{
	printf("%d", bi);
}

struct conv_arg {
	uint8_t		*buf;
	size_t		buf_sz;
	size_t		bits;

	uint8_t		m;
};

static void convolution_push_2bit(int bi, struct conv_arg *arg)
{
	size_t byte = arg->bits / 8, bit = arg->bits % 8;

	if (byte < arg->buf_sz) {
		if (bit == 0)
			arg->buf[byte] = 0;

		arg->buf[byte] |= ((bi & 0b11) << (6 - bit));
		arg->bits += 2;
	}
}

static void convolution_bit(int bi, void *private_data)
{
	struct conv_arg *arg = private_data;
	const struct fec_table *table;

	table = &state_tables[!!bi][arg->m & 0b11];

	verbose_printf(1, "bi: %d m: %d o: %d%d next_m: %d%d\n",
		       bi, arg->m,
		       bitof(table->o, 1), bitof(table->o, 0),
		       bitof(table->next_m, 1), bitof(table->next_m, 0));

	convolution_push_2bit(table->o, arg);
	arg->m = table->next_m;
}

static size_t convolution(const uint8_t *buf, size_t buf_bits, uint8_t *out_buf,
			  size_t out_buf_sz)
{
	struct conv_arg arg = {
		.buf = out_buf,
		.buf_sz = out_buf_sz,
	};

	foreach_bit(buf, buf_bits, convolution_bit, &arg);
	return arg.bits;
}

static uint8_t bufin_peek_2bit(const uint8_t *bufin, size_t *bufin_bit_idx)
{
	size_t bit_idx, byte_idx;

	bit_idx = *bufin_bit_idx;
	(*bufin_bit_idx) += 2;

	byte_idx = bit_idx / 8;
	bit_idx = bit_idx % 8;

	return (bufin[byte_idx] >> (6 - bit_idx)) & 0b11;
}

static int bufout_push_1bit(uint8_t *bufout, size_t bufout_size,
			    size_t *bufout_bit_idx, int bi)
{
	size_t bit_idx, byte_idx;

	bit_idx = *bufout_bit_idx;
	byte_idx = bit_idx / 8;
	bit_idx = bit_idx % 8;

	if (byte_idx < bufout_size) {
		if (bit_idx == 0)
			bufout[byte_idx] = 0;

		(*bufout_bit_idx) += 1;
		bufout[byte_idx] |= (bi << (7 - bit_idx));
		return 0;
	}

	return -1;
}

static uint8_t hamming_distance(uint8_t a, uint8_t b)
{
	static const uint8_t distances[] = {
		[0b00] = 0,
		[0b01] = 1,
		[0b10] = 1,
		[0b11] = 2,
	};
	uint8_t n = (a ^ b) & 0b11; /* both a and b is 2 bits number */

	return distances[n];
}

struct viterbi_trellis {
	uint8_t		bit;
	uint8_t		metric;
	uint8_t		m;
	uint8_t		prev_m;
};

struct viterbi_trellis_column {
	uint8_t				metric_min;
	struct viterbi_trellis		trellis[CONV_MAX_M];

	/* @track_back_idx: The trellis's index in this column when track back.
	 *                  setted after @viterbi_track_back.
	 *                  for debug perpose
	 */
	uint8_t				track_back_idx;
};

struct viterbi_trellis_maps {
	#define TRELLIS_COLUMN		(64 + 1) /* reserve one for head */
	size_t				n_cols;
	struct viterbi_trellis_column	cols[TRELLIS_COLUMN];
	uint8_t				bit_track_back[TRELLIS_COLUMN];
};

static void viterbi_trellis_maps_init(struct viterbi_trellis_maps *maps)
{
	memset(maps, 0xff, sizeof(*maps));

	maps->n_cols = 0;
}

static void viterbi_trellis_push(struct viterbi_trellis_column *col, uint8_t bit,
				 uint8_t distance, uint8_t m, uint8_t prev_m)
{
	struct viterbi_trellis *trellis = &col->trellis[m];

	if (distance < trellis->metric) {
		trellis->bit = bit;
		trellis->metric = distance;
		trellis->m = m;
		trellis->prev_m = prev_m;

		if (distance < col->metric_min)
			col->metric_min = distance;
	}
}

/* acs: add-compare-set */
static void viterbi_acs(struct viterbi_trellis *parent, uint8_t symbol,
			struct viterbi_trellis_column *col)
{
	for (int bit = 0; bit <= 1; bit++) {
		const struct fec_table *table = &state_tables[bit][parent->m];
		uint8_t d;

		d = hamming_distance(symbol, table->o) + parent->metric;

		viterbi_trellis_push(col, bit, d, table->next_m, parent->m);
	}
}

#define SHELL_COLOR_RESET		"\033[0m"
#define SHELL_COLOR_TEXT_RED		"\033[2;31m"
#define SHELL_COLOR_TEXT_YELLOW		"\033[2;33m"

static void print_trellis_maps(struct viterbi_trellis_maps *maps)
{
	for (int loop = 0; loop < 2; loop++) {
		printf("%s", loop == 0 ? "time: | " : "------|-");
		for (size_t col_idx = 0; col_idx < maps->n_cols; col_idx++) {
			if (loop == 0)
				printf(" %3zu  ", col_idx);
			else
				printf(" ---- ");
		}
		printf("\n");
	}

	for (int m = 0; m < CONV_MAX_M; m++) {
		printf("m %d:  | ", m);

		for (size_t col_idx = 0; col_idx < maps->n_cols; col_idx++) {
			const char *color = NULL;
			struct viterbi_trellis_column *col;
			struct viterbi_trellis *t;
			char s[32];

			col = &maps->cols[col_idx];
			t = &col->trellis[m];

			if (t->metric == 0xff) {
				printf(" ---- ");
				continue;
			}

			snprintf(s, sizeof(s), "%d/%d%c",
				 t->prev_m, t->metric, t->bit ? '.' : ' ');

			if (t->m == col->track_back_idx)
				color = SHELL_COLOR_TEXT_RED;
			else if (t->metric == col->metric_min)
				color = SHELL_COLOR_TEXT_YELLOW;

			if (color)
				printf(" %s%s%s ", color, s, SHELL_COLOR_RESET);
			else
				printf(" %s ", s);
		}
		printf("\n");
	}
}

static int viterbi_trellis_compare(const void *a, const void *b)
{
	const struct viterbi_trellis *ta = a, *tb = b;

	return ta->metric - tb->metric;
}

static int viterbi_track_back(struct viterbi_trellis_maps *maps, uint8_t *last_m)
{
	struct viterbi_trellis_column ep; /* sorted endpoint */
	struct viterbi_trellis *t;
	int col_idx, bit_idx = 0;

	/* sort it */
	memcpy(&ep, &maps->cols[maps->n_cols - 1], sizeof(ep));
	qsort(ep.trellis, CONV_MAX_M, sizeof(ep.trellis[0]),
	      viterbi_trellis_compare);

	/* ep.trellis[0] is the best one */
	t = &ep.trellis[0];
	col_idx = maps->n_cols - 1;
	memset(maps->bit_track_back, 0, sizeof(maps->bit_track_back));

	if (last_m)
		*last_m = t->m;

	while (col_idx > 0) {
		uint8_t prev_m = t->prev_m;

		maps->bit_track_back[bit_idx++] = t->bit;
		maps->cols[col_idx].track_back_idx = t->m;

		--col_idx;
		t = &maps->cols[col_idx].trellis[prev_m];
	}

	if (option_verbose > 1) {
		printf("Viterbi trackback path:\n");
		print_trellis_maps(maps);
	}

	return 0;
}

struct viterbi_decode_arg {
	uint8_t			m;

	const uint8_t		*in;
	size_t			in_bits;
	size_t			in_bit_idx;

	uint8_t			*out;
	size_t			out_buf_sz;
	size_t			out_bit_idx;
};

static ssize_t viterbi_decode_part(struct viterbi_decode_arg *arg)
{
	struct viterbi_trellis_maps maps;

	viterbi_trellis_maps_init(&maps);

	/* push the first entry to maps. */
	memset(&maps.cols[0].trellis[0], 0, sizeof(struct viterbi_trellis));
	maps.cols[0].trellis[0].m = arg->m & 0b11;
	maps.cols[0].metric_min = 0;
	maps.n_cols = 1;

	for (; maps.n_cols < TRELLIS_COLUMN && arg->in_bit_idx < arg->in_bits;
	     maps.n_cols++) {
		struct viterbi_trellis_column *col = &maps.cols[maps.n_cols];
		uint8_t symbol = bufin_peek_2bit(arg->in, &arg->in_bit_idx);

		for (size_t j = 0; j < CONV_MAX_M; j++) {
			/* the previous column and in the same line */
			struct viterbi_trellis *prev = &col[-1].trellis[j];

			if (prev->metric == 0xff)
				continue;

			viterbi_acs(prev, symbol, col);
		}
	}

	if (option_verbose) {
		printf("Viterbi trellis maps:\n");
		print_trellis_maps(&maps);
	}

	/* track back and push the last_m to arg */
	viterbi_track_back(&maps, &arg->m);

	for (int n = maps.n_cols - 1 - 1; n >= 0; n--)
		bufout_push_1bit(arg->out, arg->out_buf_sz,
				 &arg->out_bit_idx,
				 maps.bit_track_back[n]);

	return arg->out_bit_idx;
}

static ssize_t viterbi_decode(const uint8_t *buf, size_t buf_bits,
			      uint8_t *out_buf, size_t out_buf_sz)
{
	struct viterbi_decode_arg arg = {
		.m = 0,
		.in = buf,
		.in_bits = buf_bits / 2 * 2, /* 1bit generate 2bits conv code */
		.out = out_buf,
		.out_buf_sz = out_buf_sz,
	};

	while (arg.in_bit_idx < arg.in_bits) {
		ssize_t ret;

		ret = viterbi_decode_part(&arg);
		if (ret < 0)
			return ret;
	}

	return arg.out_bit_idx;
}

enum {
	OPTION_DECODE,
	OPTION_ENCODE,
	OPTION_VERBOSE,
};

static struct option long_options[] = {
	/* name		has_arg,		*flag,	val */
	{ "decode",	no_argument,		NULL,	OPTION_DECODE	},
	{ "encode",	no_argument,		NULL,	OPTION_ENCODE	},
	{ "verbose",	no_argument,		NULL,	OPTION_VERBOSE	},
	{ NULL,		0,			NULL,	0		},
};

static void print_usage(void)
{
	fprintf(stderr, "Usage: viterbi [OPTIONS] \"binary-strings\"\n");
	fprintf(stderr, "-d --decode:            decode it\n");
	fprintf(stderr, "-e --encode:            encode it\n");
}

int main(int argc, char **argv)
{
	uint8_t input[128] = { 0 }, output[sizeof(input) * 2] = { 0 };
	size_t input_bits = 0;
	const char *endp;
	int ret, decode = 0;

	while (1) {
		int option_index = 0;
		int c;

		c = getopt_long(argc, argv, "vde", long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
		case OPTION_DECODE:
		case OPTION_ENCODE:
			decode = c == OPTION_DECODE;
			break;
		case 'd':
		case 'e':
			decode = c == 'd';
			break;
		case 'v':
		case OPTION_VERBOSE:
			option_verbose++;
			break;
		}
	}

	if (!(optind < argc)) {
		print_usage();
		return -1;
	}

	input_bits = str01_to_buffer(argv[optind], &endp, input,
				     sizeof(input), 0);
	if (input_bits == 0 || *endp != '\0') {
		fprintf(stderr, "Wrong binary stream after <%s>\n", endp);
		return -1;
	}

	if (decode)
		ret = viterbi_decode(input, input_bits, output, sizeof(output));
	else
		ret = convolution(input, input_bits, output, sizeof(output));

	if (ret < 0)
		return ret;

	if (ret > 0) {
		foreach_bit(output, ret, print_bit, NULL);
		printf("\n");
	}

	return 0;
}
