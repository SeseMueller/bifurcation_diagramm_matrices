# Calculating and painting the Bifurcation Diagram

## (But not the way you think you know it)

I have over the years written a few programs with the goal of visualizing the bifurcation diagram of the logistic map.
But between switching language, doing microoptimizations and parallelizing the code, I have never
really written a program that I was happy with. They are all slow eventually or become
inaccurate at some point. So I wanted to know whether it is possible to calculate the logistic map in
sub - `O(ITER)` time.

I am sad to say that this program is no different (Slow and inaccurate, that is).
_BUT_, it uses two different approaches than the convential one, and I think that is worth
something.

# The usual approach

The usual approach to visualizing the logistic map is to fix the `r` parameter, iterate the map a bit
and then record the next `n` iterations. This is repeated for all desired `r` values. The results
are then plotted in a scatter plot. This has the advantage of being very simple to implement and
quite fast, scaling linearly in both iterations and output picture size. There are only a few optimizations
that can be done to speed it up without changing the algorithm. One would be to go closer to the
metal, using SIMD instructions, parallelizing the code, etc.
However, this approach often needs an increadibly large amount of iterations to get a good picture.

# The first approach - buckets: O(ITERS)

My first idea was to use "buckets", inspired by the [bucket sort](https://en.wikipedia.org/wiki/Bucket_sort) algorithm.
The idea is to divide the `[0, 1]` interval into `n` buckets, and then not map single values to the buckets,
but rather intervals. This way, we can iterate over the buckets and calculate the next iteration for all
values in the bucket at once.
To implement this, I calculate a matrix that, for a fixed `r` value, tells me the fraction of values that will
go from one bucket to another per iteration. Then, for a starting value, I can calculate the fraction of values
that will end up in each bucket after `n` iterations. This is done by multiplying the bucket vector with the
matrix `n` times. This has the amazing property that there isn't a single number that is being iterated over,
but rather an "infinite amount" (but not really) of numbers. This means that we can get a good picture with
a relatively small amount of iterations. The downside is that we need to calculate the matrix, as well as
multiply it with the bucket vector `n` times. This is done in `O(BUCKETS^3)` time.
By the way, if you want render the image, _every pixel has to be a bucket_ or you lose resolution.

This is already slower than the usual approach, but can we, given my original goal of calculating the logistic
map in sub - `O(ITER)` time, do better?

The Answer is yes.

# The second approach - the matrix: O(log(ITER))

The second approach is very similar to the first, but instead of multiplying the matrix with the bucket vector
`n` times, we can multiply the matrix with **itself** `log(ITER)` times. This is possible because the matrix
is a [stochastic matrix](https://en.wikipedia.org/wiki/Stochastic_matrix), meaning that all its entries are
non-negative and the sum of each column is 1. (Also because we did it before, whether we calculate ((A\*A)^4)\*x
or x <= Ax 16 times doesn't make a difference mathematically).

**This way, we can actually calculate the logistic map in sub - `O(ITER)` time! My dream come true!**

But, as you might have noticed, multiplying a matrix with itself is increadibly slow (I think the underlaying implementation is O(n^3)). This means that the second approach is actually slower than the first one. But, it is still a cool idea.

### Getting side tracked

At this point I looked into a few other things. Maybe it's possible (and hopefully faster) to calculate the matrices
beforehand? The answer is that this works and is marginally faster, but also takes up a lot of memory (I don't have 64 GB memory).

Another Idea I had was to do fewer iteration when the values don't change much. This can easily be done by finding the
maximum value in the vector and breaking once it reaches below a certain threshold. This is a good idea and makes the
program a bit faster.

Higher thresholds can lead to faster times, but also to less accurate pictures. I found that a threshold of 0.001 is a
good balance.

But enough with the side tracking, let's get back to the main topic.

Is it possible to faster? (In terms of O-time relative to iterations)

The answer is yes, but it's not pretty.

# The third approach - heresy: O(1)

I was writing the first few approaches, sleep deprived and with appearently not enough knowledge about linear algebra
and markov chains. I was thinking about what I was actually doing and tried to rearrange the equation.

With our Matrix $A$ and our vector $x$, we are calculating $A^n*x$. But we want $Ax=x$. Using linear algebra, we can
rearrange:

$$
Ax=x \\
Ax-x=0 \\
Ax-Ix=0 \\
(A-I)x=0 \\
x=0
$$

As it turns out, the "stable point" that is reached stochastically by iterating the logistic map ad infinitum in the
usual approach by the law of large numbers actually doesn't exist here.

BUT, it just doesn't exist mathematically. As suckerpinch pointed out [in his amazing video](https://www.youtube.com/watch?v=Ae9EKCyI1xU), the floating point errors slightly violate the conventions of mathematics. This means that, instead of searching for an $x$ such that $Ax=x$, we can search for an $x$ such that $Ax=x+\epsilon$ for some small $\epsilon$. In the previous approaches, we did this by iterating the matrix and hoping it will converge. As we have seen, that is actually the case, so we must be able to find such an almost-stable point.

So I had the totally brilliant idea to just do that:

$$
Ax=x+\epsilon \\
Ax-x=\epsilon \\
(A-I)x=\epsilon \\
x=\epsilon(A-I)^{-1}
$$

Using the fact that the sum of all entries in $x$ is 1, we can actually calculate $x$ exactly. This is done by
calculating the inverse of $A-I$ and multiplying it with $\epsilon$. This is done in `O(BUCKETS^3)` time, but that doesn't bother us anymore. We can now calculate the logistic map in `O(1)` time, not even needing to specify iteration count!

But as you might have realized, there is a small, tiny, little problem (Besides having a quite slow program). We need to calculate the inverse of $A-I$.

Which doesn't have to exist.

## Heresy 2 - Electric Boogaloo: O(1)?

If you look at the problem at hand, you might realize that we have a bit more freedom to the matrix than we thought.
Because we can technically treat the values in the buckets as points, we can let them "stay" where they are. This
means that we can (as lonng as we preserve the fact that the sum of all entries in a column is 1) increase the values
on the diagonal of the matrix _by any amount we want_. This means that we can just _make_ the matrix invertible.

The way I chose to do it was to multiply the entire matrix with $0.8$ and then add $0.2$ to the diagonal. This way,
most matrices become invertible! We have finally found a way to reliably calculate the logistic map in `O(1)` time!

But have we really?

## Heresy 3 - The iteration strikes back: O($\infty$)???

Yes, only most matrices become invertible. Some, by chance, don't. We _could_ change the matrix diagonal by slight amounts, but that would still not work definitively. (Also I just didn't think of it at the time)

The solution is obvious: We just pertubate the matrix a bit. What does it matter if we change the matrix by a tiny amount? It's still the same matrix, right?

Surprisingly, adding a small amount of noise (<0.0005) to the matrix actually makes it invertible. And if it doesn't,
we just add a bit more noise.
The only downside is that the resulting vector can, if it was close to 0 before, become negative. This is easily fixed by setting all negative values to 0 and then normalizing the vector.

But the price we pay for this is that we have sacrificed our `O(1)` worst-case time, because the matrix might just
never become invertible. This is however very unlikely. In fact, the randomness gets called less than 1 in 2000 matrix pseudo-inversions.

# The final result

In this repository, you can find my three approaches to calculating the logistic map. The code is not very good, but it
works. The first approach is obviously the fastest, but the others also work. The memory usage is not very high.

Because I tried to optimize the code a bit for my hardware, it might not work for you out of the box. If you, for some
reason, decide to run this locally, you might want to change the Feature layout in Cargo.toml. Just toggle the
commented and uncommented lines.

# Conclusion

The main takeaway from this is that Big-O notation is not everything, especially if you are using multiple variables in
the equation. The different approaches have different Big-O times, both in terms of iterations and buckets.

Another thing is that you really do need a mathematical background to do this kind of stuff. I was lucky enough to have
had a few years of linear algebra in uni, but I should brush up on my markov chains. There is probably a better way to find a stable point, since you could represent the matrix as a directed weighted graph and then use some markov chain algorithm to find the stable point.

Lastly, you should think about what different algorithms you think of actually do, i.E. what the computer is doing when
it executes your code. I thought that eliminating the need to iterate a lot would be a good idea and that computers are
optimized to do matrix multiplication. But I didn't think about the fact that the matrix is actually quite big and that
every pixel needs to be a bucket.

Nontheless, I am happy to have been wrong, because I learned quite a bit about applied linear algebra and high performance rust.
