---
title: Test Quarto with Python
date: 2024-06-09T00:00:00.000Z
author: Sean Hackett
layout: post
tags:
  - test
  - python
  - quarto
jupyter: python3
---


# Quarto Python Test

import math import random

# Simple calculations

numbers = \[1, 2, 3, 4, 5, 6, 7, 8, 9, 10\] squares = \[x\*\*2 for x in
numbers\]

print(“Numbers:”, numbers) print(“Squares:”, squares) print(“Sum of
squares:”, sum(squares))

# Generate some random data

random.seed(42) random_data = \[random.randint(1, 100) for \_ in
range(10)\] print(“Random data:”, random_data) print(“Average:”,
sum(random_data) / len(random_data))
