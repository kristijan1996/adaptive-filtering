# adaptive-filtering
Implementation of adaptive filter class in Matlab/Octave

## Intro

In real systems there is often a need for processing signals which have unknown characteristics or are variable in time. In these systems it is practically impossible to use standard filters and digital signal processing algorithms, because their parameters have to be constantly up-to-date to match the characteristics of the concrete problem.
In this project, a class of adaptive filters is implemented, which could provide useful tool for filtering signals and its usage will be explained in the following sections.

## The main idea

While designing any general time-varying (adaptive) system, the idea is to adjust the system (filter) coefficients, i.e. system's transfer functions, iteratively so that the system produces the output closest to the desired one. Adaptive filters, as well as neural networks in their basic form, do this by minimizing the error between the desired output and the current filter output through iterations.
On the following image a block diagram of adaptive filter is depicted.

<p align="center">
  <img width="460" src="/images/1.png">
</p>

_x(n)_ is the input signal, _y(n)_ output obtained by filtering _x_, _d(n)_ the desired filter output and _e(n)_ the difference between the desired and current output. The training algorithm adjusts the filter coefficients in order to minimize the mean squared error.

## Design

Typically, adaptive filters are _FIR_ filters, because of the demand for stability and linear phase response. On the following image, a block diagram of _FIR_ filter is shown.

<p align="center">
  <img width="460" src="/images/2.png">
</p>

Basically, this structure allows for estimating the output sample _y(n)_ by processing input samples _x(n), x(n-1), ... x(n - L+1)_, where _L_ is the order of the filter. This estimation of output by observing current and previous inputs to the systems renders it causal, which means it can work in real-time if implemented, let's say, in _C_ on a microcontroller.
The depicted _FIR_ structure boils down to simple matrix (vector) multiplication in Matlab.

## Use cases

One of the many potential use cases of adaptive filters is in estimation of an unknown system's transfer function, as shown on the following diagram.

<p align="center">
  <img width="460" src="/images/systemIdentification.png">
</p>

By setting _d(n)_ to be the output of the unknown system, adaptive filter tunes its coefficients through iterations in order to resemble _P(z)_ from input-output point of view, regardless of the inner structure of _P(z)_. 

## Example

Let's define a pure _10 Hz_ sinusoid as our desired output from an adaptive filter which is fed with a noisy sine wave additionally polluted with unwanted sinusoid interference on _50 Hz_. 

```
fs = 1000;
dt = 1/fs;
f = 10;
t = 0:dt:5-dt;
x_d = sin(2*pi*f*t);
x_in = x_d + sin(2*pi*(f+40)*t) + 0.5*randn(1, numel(x_d));
```
By instantiating the object of AdaptFilter class of order _100_ and learning rate _0.01_, we can call _LMS_ method on it, by providing it with noisy input _x\_in_ and desired ouput _x\_d_.

```
filter = AdaptFilter(100, 0.01);

[y, e, h] = filter.LMS(x_in, x_d);
```

On the following graphs we can see the noisy input as well as ouputs from the filter we implemented and the adaptive filter built into Matlab.

<p align="center">
  <img width="600" src="/images/filterOut.png">
</p>

The noisy sine wave has been cleared of interferences, which can also be seen from signal spectrums.

<p align="center">
  <img width="600" src="/images/outSpectrumIn.png">
</p>
<p align="center">
  <img width="600" src="/images/outSpectrumClean.png">
</p>

# Author
Kristijan Mitrovic

# License
Free in every way
