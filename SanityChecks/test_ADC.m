 clear all

[data, timestamps, info] = load_open_ephys_data_faster('100_ADC1.continuous');
subplot(3,3,1)
plot(timestamps,data)

[data2, timestamps2, info2] = load_open_ephys_data_faster('100_ADC2.continuous');
subplot(3,3,2)
plot(timestamps2,data2)
%ylim([-4.5 -4])

[data3, timestamps3, info3] = load_open_ephys_data_faster('100_ADC3.continuous');
subplot(3,3,3)
plot(timestamps3,data3)
%ylim([-4.5 -4])

[data4, timestamps4, info4] = load_open_ephys_data_faster('100_ADC4.continuous');
subplot(3,3,4)
plot(timestamps4,data4)
ylim([-4.5 -4])

[data5, timestamps5, info5] = load_open_ephys_data_faster('100_ADC5.continuous');
subplot(3,3,5)
plot(timestamps5,data5)
ylim([-4.5 -4])

[data6, timestamps6, info6] = load_open_ephys_data_faster('100_ADC6.continuous');
subplot(3,3,6)
plot(timestamps6,data6)
ylim([-4.5 -4])

[data7, timestamps7, info7] = load_open_ephys_data_faster('100_ADC7.continuous');
subplot(3,3,7)
plot(timestamps7,data7)
ylim([-4.5 -4])

[data8, timestamps8, info8] = load_open_ephys_data_faster('100_ADC8.continuous');
subplot(3,3,8)
plot(timestamps8,data8)
ylim([-4.5 -4])
