function max_=row_max(a)
    max_ = -inf;
    [m,n] = size(a);
    for i=1:m
        now_row=a(i,:);
        now_sum=sum(now_row);
        disp(now_sum);
        if now_sum>max_
            max_=now_sum;
        end
    end
end
