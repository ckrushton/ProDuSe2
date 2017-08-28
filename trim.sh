#!/usr/bin/env bash
./trim.py --config:batch: "samples[].[merge({output: join(\`\`, [output, \`trimmed/\`, name, \`.fastq.gz\`]), input:fastqs[0], mate: fastqs[1]}, trim)] | []" "$@"