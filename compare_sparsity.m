filename = 'yeast';
str1 = [filename '.mtx'];
str2 = ['output_' filename '.txt'];

original = load(str1);
original = original(2:end,:);
new = load(str2);

S = spconvert(original);
S2 = spconvert(new);

subplot(1,2,1); spy(S);
subplot(1,2,2); spy(S2);

