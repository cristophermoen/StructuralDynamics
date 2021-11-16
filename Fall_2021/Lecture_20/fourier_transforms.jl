
#' # Lecture 20 - Fourier transforms

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

#' It is sometimes interesting and useful to calculate the frequency content of a signal.  


using FFTW

#' Followed this example.
#' https://medium.com/@khairulomar/deconstructing-time-series-using-fourier-transform-e52dd535a44e

#' Define a sine wave.
f = 10 #cycles per second, Hz
ω = 2 * π * f
Fs = 1000  #samples per second , dt = 0.001 seconds
t = range(0.0, 1.0, step = 1/Fs)
p = 1.0 * sin.(ω * t)

plot(t, p)

signal_fft = fft(p)

n = length(t)

frequency = Fs/2 .* range(0.0, 1.0, length = floor(Int, n/2))

y_m = 2/n * abs.(signal_fft[1:length(frequency)])

plot(frequency, y_m, xlims = (0.0, 20.0))


#' Now let's try this same thing with an earthquake ground motion.

#' Let's consider an earthquake ground motion from an earthquake in 1979 centered in Imperial Valley, California.  The ground accelerations are obtained from the [Center for Engineering Strong Motion Data](https://www.strongmotioncenter.org/cgi-bin/CESMD/iqr_dist_DM2.pl?IQRID=ImperialValley79&SFlag=0&Flag=3).  

#' Read in ground acceleration data, note dims=(rows,columns) in text file that has 711 lines and 8 columns.
using DelimitedFiles
earthquake_data = readdlm("/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Lecture_12/Imperial_Valley_1979.txt",Float64,dims=(711,8))

#' Rearrange ground accleration data as a vector.
earthquake_data=transpose(earthquake_data)   #transpose matrix
utt_g=reshape(earthquake_data,(5688,1))  #reshape matrix into a vector, note 8*711=5688
utt_g=vec(utt_g)   #tell Julia to make the data a 1D vector

#' Units come in as cm/sec^2 from strongmotion.org, let's change them to m/sec^2.
utt_g=utt_g ./ 100.0

#' Define time range of the earthquake assuming the first acceleration reading occurs at t=0.0 seconds.  The ground acceleration is provided every 0.010 seconds.
t_eq=collect(range(0,length=5688,stop=5687*0.01))   #total time is 5687*0.01 seconds

#' Plot ground acceleration over time.
using Plots
plot(t_eq,utt_g,linewidth=1,title="Imperial Valley 1979 earthquake",
    xaxis="time, t (sec.)",yaxis="ground acc., m/sec^2",legend=false)


Fs = 1/0.01    
eq_signal_fft = fft(utt_g)
n = length(t_eq)
frequency = Fs/2 .* range(0.0, maximum(t_eq), length = floor(Int, n/2))
fourier_coeffs = 2/n * abs.(eq_signal_fft[1:length(frequency)])
plot(frequency, fourier_coeffs, xlims=(0.0, 100.0))