
invested_money = 7500;
for n = 1:20
    invested_money = invest_money(invested_money);
end

disp(invested_money)

function new_money = invest_money(old_money)

    new_money = 7500+1.08*(old_money);
end
