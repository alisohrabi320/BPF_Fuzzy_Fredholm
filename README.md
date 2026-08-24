# BPF_Fuzzy_Fredholm
ساختار مخزن (Repository)
LIFFIE-BPF-Solver/
│
├── README.md                       ← توضیح فارسی پروژه
├── requirements.txt                ← پیش‌نیازهای پایتون
├── LICENSE                         ← مجوز (اختیاری)
│
├── code/
│   ├── full_solver.py              ← کد کامل (۶ سکشن)
│   ├── section1_base.py            ← سکشن ۱: توابع پایه
│   ├── section2_examples.py        ← سکشن ۲: حل‌کننده مثال‌ها
│   ├── section3_errors.py          ← سکشن ۳: محاسبه خطاها
│   ├── section4_convergence.py     ← سکشن ۴: نرخ همگرایی
│   ├── section5_analysis.py        ← سکشن ۵: مقایسه و کنترل
│   └── section6_plots.py           ← سکشن ۶: نمودار و main
│
├── outputs/
│   ├── example1_tables.xlsx
│   ├── example2_tables.xlsx
│   ├── example3_tables.xlsx
│   ├── all_examples_convergence.png
│   ├── all_examples_convergence.pdf
│   ├── center_cell_convergence.png
│   ├── center_cell_convergence.pdf
│   └── final_report.txt
│
└── docs/
    ├── usage_fa.md                 ← راهنمای فارسی
    └── algorithm_description.md    ← توضیح الگوریتم‌ها
# حل عددی معادلات انتگرال فردهلم فازی شهودی با توابع بلوک-پالس

این مخزن شامل کدهای پایتون پایان‌نامه کارشناسی ارشد با عنوان «حل عددی معادلات انتگرال فردهلم فازی شهودی با استفاده از توابع بلوک-پالس و تعمیم به حالت غیرخطی» است.

---

## 📋 ویژگی‌ها

- ✅ حل سه مثال عددی (دو خطی + یک غیرخطی)
- ✅ پیاده‌سازی روش BPF با قانون ضرب فازی
- ✅ حل دستگاه غیرخطی با تکرار نقطه ثابت
- ✅ محاسبه خطا و نرخ همگرایی
- ✅ ذخیره جداول در Excel
- ✅ رسم و ذخیره نمودارهای همگرایی (PNG و PDF)

---

## 🧩 مثال‌های پیاده‌سازی‌شده

| مثال | نوع | هسته | غیرخطی؟ |
|:---:|:---:|:---:|:---:|
| ۱ | خطی | \(k(t,s)=ts\) | ❌ |
| ۲ | خطی | \(k(t,s)=t^2s\) | ❌ |
| ۳ | غیرخطی | \(k(t,s)=ts\) | ✅ \(g(y)=y^2\) |

---

## ⚙️ نصب پیش‌نیازها

`bash
pip install numpy pandas scipy matplotlib openpyxl
🚀 اجرا

اجرای کامل (هر ۴ سکشن در یک فایل):

bash
python code/full_solver.py


اجرای مثال خاص:

برنامه از شما می‌پرسد کدام مثال را حل کنید:


کدام مثال را حل کنم؟
  ۱: مثال ۱ (خطی، k=ts)
  ۲: مثال ۲ (خطی، k=t²s)
  ۳: مثال ۳ (غیرخطی، y²)
  ۴: هر سه
انتخاب (۱-۴):


---

📂 خروجی‌ها

پس از اجرا، خروجی‌های زیر تولید می‌شوند:

فایل توضیح
example{1,2,3}_tables.xlsx جداول خطا برای m=8 و m=16
example{1,2,3}_convergence.png نمودار همگرایی هر مثال
example{1,2,3}_convergence.pdf نسخه PDF نمودار
all_examples_convergence.png نمودار ترکیبی هر سه مثال
all_examples_convergence.pdf نسخه PDF ترکیبی

---

🧮 ساختار کد 

سکشن فایل توضیح
full_solver.py
│
├── سکشن ۱: ایمپورت و توابع پایه BPF
│   ├── import numpy, pandas, scipy, matplotlib
│   ├── kernel_ts()          → هسته k(t,s) = ts
│   ├── kernel_t2s()         → هسته k(t,s) = t²s
│   ├── compute_K_matrix()   → ساخت ماتریس K با انتگرال دوگانه
│   ├── solve_linear()       → حل دستگاه خطی + قانون ضرب فازی
│   └── reconstruct()        → بازسازی جواب از ماتریس Y
│
├── سکشن ۲: حل‌کننده‌های سه مثال
│   ├── solve_example1()     → مثال ۱ (خطی، k=ts)
│   ├── solve_example2()     → مثال ۲ (خطی، k=t²s)
│   ├── g()                  → تابع غیرخطی g(y) = y²
│   └── solve_example3()     → مثال ۳ (غیرخطی، تکرار نقطه ثابت)
│
├── سکشن ۳: محاسبه خطاها
│   ├── make_table()         → ساخت جدول خطا برای هر t
│   ├── compute_max_error()  → خطای Full-Grid
│   └── compute_center_error() → خطای Center-Cell
│
├── سکشن ۴: همگرایی و نرخ‌ها
│   ├── run_example()        → اجرای یک مثال و جمع‌آوری خطاها
│   ├── compute_rate()       → نرخ همگرایی Full-Grid
│   └── compute_center_convergence() → نرخ Center-Cell
│
├── سکشن ۵: جدول مقایسه و کنترل خودکار
│   ├── analyze_all_errors() → مقایسه Full-Grid و Center-Cell
│   └── automatic_checks()   → بررسی کاهش خطا و منطقی بودن نرخ
│
└── سکشن ۶: رسم نمودارها و اجرا
    ├── plot_convergence_all() → نمودار ترکیبی log-log
    ├── plot_center_convergence() → نمودار Center-Cell
    ├── save_report()          → ذخیره گزارش نهایی
    └── main()                 → اجرای همه مراحل

---

📊 روش محاسبه خطا

خطای بیشینه با فرمول زیر محاسبه می‌شود:

\text{Max Error} = \max |y_{\text{approx}} - y_{\text{exact}}|

روی ۱۸ نقطه (t = 0, 0.125 و r = 0, 0.125, \ldots, 1) و برای هر ۴ مؤلفه.

---

🎯 نرخ همگرایی

نرخ همگرایی از رابطه زیر به دست می‌آید:

\text{Rate} = \frac{\log(e_1/e_2)}{\log(m_2/m_1)}

---

📚 منابع

· Ghazanfari, B. (2025). NIFS, 31(1), 99–110.
· Zadeh, L. A. (1965). Information and Control, 8, 338–353.
· Atanassov, K. (1986). Fuzzy Sets and Systems, 20(1), 87–96.

---

📝 مجوز

این پروژه برای اهداف آموزشی و پژوهشی منتشر شده است.



---

# requirements.txt



numpy>=1.21
pandas>=1.3
scipy>=1.7
matplotlib>=3.4
openpyxl>=3.0



---

# docs/usage_fa.md


markdown
# راهنمای استفاده

## گام ۱: نصب

در ترمینال:

bash
pip install numpy pandas scipy matplotlib openpyxl


گام ۲: کپی کد

هر ۴ سکشن را به ترتیب در یک فایل full_solver.py کپی کنید.

گام ۳: اجرا

bash
python full_solver.py


گام ۴: انتخاب مثال

برنامه از شما می‌پرسد کدام مثال حل شود. عدد ۱ تا ۴ را وارد کنید.

گام ۵: دریافت خروجی

فایل‌های Excel و PNG/PDF در همان پوشه ذخیره می‌شوند.
