**Results Summary:**

We have chosen to implement (1) a PID controller, (2) an LQR controller, and (3) a Luenberger Observer. 

We report the trends for our simulation results below:

![scoreVsAmpl_sine](https://github.com/user-attachments/assets/8ee3f203-004d-4170-8039-61973e94c6c0)

![scoreVsPeriod_square](https://github.com/user-attachments/assets/59e67a05-1710-4332-bf57-1dc0de7ea1ef)

![scoreVsAmpl_square](https://github.com/user-attachments/assets/b4eed003-09c5-482c-8661-fbde873f97ea)

![scoreVsPeriod_square](https://github.com/user-attachments/assets/d3be6811-ccaa-4e37-a629-558372ce1ba6)

To better understand how each controller performs, we compare scores across various amplitudes and periods. First, we hold the period constant. For the sine wave, we see a nearly linear increase in score for LQR and a slightly worse than linear increase for PID. For the square wave, we see a more exponential increase, where the score starts to blow up for both LQR and PID. Next, we hold the amplitude constant. For the sine wave, we see a nearly linear decrease in score for PID and a slightly faster decrease for LQR. For the square wave, we see a fairly flat decline in score. However, the PID has an overall higher scores. Overall, we see slightly better score trends for the LQR controller. 


Further, we give a brief comparative analysis when comparing single runs:

The highest magnitude of the control input is greater for PID compared to LQR. 
The PID performs slightly better than LQR in terms of the score for both sine and square waves. 
Both PID and LQR tend to do better for sine waves compared to square waves in terms of both tracking performance and energy cost. 

Screen recording of the controller animation:

**LQR**

https://github.com/user-attachments/assets/d2121167-b21d-4345-910c-66ce1deb9571


https://github.com/user-attachments/assets/705f7859-934f-4d94-ade3-8f627197541a



**PID**

https://github.com/user-attachments/assets/9bfe72c4-3f3e-4e41-bc4e-a8fecc587da3


https://github.com/user-attachments/assets/5c26c1ae-1474-46c9-917f-32cfcda8475b


