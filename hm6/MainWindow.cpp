#include "MainWindow.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QPainter>
#include <QRegularExpression>
#include <QMessageBox>

PlotWidget::PlotWidget(QWidget *parent) : QWidget(parent) {
    setMinimumSize(500, 500);
}

void PlotWidget::setGenerators(const QVector<Monom>& gens) {
    m_generators = gens;
    update();
}

void PlotWidget::setMode(bool showIdeal) {
    m_showIdeal = showIdeal;
    update();
}

void PlotWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event);
    
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    
    int w = width();
    int h = height();
    int marginLeft = 50;
    int marginBottom = 40;
    int cellW = (w - marginLeft - 20) / (MAX_COORD + 1);
    int cellH = (h - marginBottom - 20) / (MAX_COORD + 1);
    
    int startX = marginLeft;
    int startY = h - marginBottom;
    
    painter.setPen(QPen(Qt::black, 2));
    painter.drawLine(startX - 5, startY, w - 20, startY);
    painter.drawLine(startX, startY + 5, startX, 20);
    
    painter.drawText(w - 30, startY - 5, "m");
    painter.drawText(startX + 5, 25, "n");
    
    for (int m = 0; m <= MAX_COORD; m++) {
        for (int n = 0; n <= MAX_COORD; n++) {
            int x = startX + m * cellW;
            int y = startY - (n + 1) * cellH;
            
            painter.setPen(QPen(Qt::lightGray, 0.5));
            painter.drawRect(x, y, cellW, cellH);
            
            bool inIdeal = false;
            for (const Monom& gen : m_generators) {
                if (m >= gen.x && n >= gen.y) {
                    inIdeal = true;
                    break;
                }
            }
            
            if (m_showIdeal) {
                painter.fillRect(x + 1, y + 1, cellW - 1, cellH - 1,
                                 inIdeal ? Qt::green : Qt::red);
            } else {
                painter.fillRect(x + 1, y + 1, cellW - 1, cellH - 1,
                                 !inIdeal ? Qt::green : Qt::red);
            }
        }
    }
    
    painter.setPen(QPen(Qt::black, 1));
    for (int m = 0; m <= MAX_COORD; m++) {
        int x = startX + m * cellW + cellW / 2 - 5;
        int y = startY + 15;
        painter.drawText(x, y, QString::number(m));
    }
    
    for (int n = 0; n <= MAX_COORD; n++) {
        int x = startX - 25;
        int y = startY - (n + 1) * cellH + cellH / 2 + 5;
        painter.drawText(x, y, QString::number(n));
    }
}

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    QWidget* centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    
    QVBoxLayout* mainLayout = new QVBoxLayout(centralWidget);
    
    QGroupBox* inputGroup = new QGroupBox("Введите образующие мономы");
    QVBoxLayout* inputLayout = new QVBoxLayout(inputGroup);
    
    QLabel* hintLabel = new QLabel(
        "Пример: x^6, x^2y^3, xy^7\n"
        "Формат: x^a y^b (порядок букв не важен)");
    hintLabel->setWordWrap(true);
    inputLayout->addWidget(hintLabel);
    
    m_inputEdit = new QLineEdit();
    m_inputEdit->setText("x^6, x^2y^3, xy^7");
    inputLayout->addWidget(m_inputEdit);
    
    m_generateButton = new QPushButton("Построить");
    inputLayout->addWidget(m_generateButton);
    
    mainLayout->addWidget(inputGroup);
    
    QGroupBox* modeGroup = new QGroupBox("Режим отображения");
    QHBoxLayout* modeLayout = new QHBoxLayout(modeGroup);
    
    m_modeCombo = new QComboBox();
    m_modeCombo->addItem("(a) Мономы, принадлежащие идеалу");
    m_modeCombo->addItem("(b) Мономы, которые могут быть в остатке");
    modeLayout->addWidget(m_modeCombo);
    
    mainLayout->addWidget(modeGroup);
    
    m_plotWidget = new PlotWidget();
    mainLayout->addWidget(m_plotWidget);
    
    m_statusLabel = new QLabel("Введите образующие и нажмите 'Построить'");
    mainLayout->addWidget(m_statusLabel);
    
    connect(m_generateButton, &QPushButton::clicked, this, &MainWindow::onGenerateClicked);
    connect(m_modeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &MainWindow::onModeChanged);
    
    onGenerateClicked();
}

QVector<Monom> MainWindow::parseGenerators(const QString& text) {
    QVector<Monom> result;
    QString clean = text;
    clean.remove(' ');
    clean.remove('*');
    
    QStringList parts = clean.split(',', Qt::SkipEmptyParts);
    
    for (const QString& part : parts) {
        Monom mon = {0, 0};
        
        QRegularExpression reX("x\\^(\\d+)");
        QRegularExpression reY("y\\^(\\d+)");
        QRegularExpression reXsimple("x(?!\\^)");
        QRegularExpression reYsimple("y(?!\\^)");
        
        QRegularExpressionMatch matchX = reX.match(part);
        QRegularExpressionMatch matchY = reY.match(part);
        
        if (matchX.hasMatch()) {
            mon.x = matchX.captured(1).toInt();
        } else if (reXsimple.match(part).hasMatch()) {
            mon.x = 1;
        }
        
        if (matchY.hasMatch()) {
            mon.y = matchY.captured(1).toInt();
        } else if (reYsimple.match(part).hasMatch()) {
            mon.y = 1;
        }
        
        if (mon.x == 0 && mon.y == 0) {
            bool ok;
            part.toInt(&ok);
            if (!ok && !part.isEmpty()) {
                QMessageBox::warning(this, "Ошибка", "Не удалось разобрать моном: " + part);
            }
            continue;
        }
        
        result.append(mon);
    }
    
    return result;
}

void MainWindow::onGenerateClicked() {
    QString input = m_inputEdit->text();
    QVector<Monom> gens = parseGenerators(input);
    
    if (gens.isEmpty()) {
        m_statusLabel->setText("Ошибка: не удалось распознать ни одного образующего монома");
        return;
    }
    
    QString gensStr;
    for (const Monom& mon : gens) {
        gensStr += QString("x^%1 y^%2, ").arg(mon.x).arg(mon.y);
    }
    gensStr.chop(2);
    m_statusLabel->setText("Образующие: " + gensStr);
    
    m_plotWidget->setGenerators(gens);
}

void MainWindow::onModeChanged(int index) {
    m_plotWidget->setMode(index == 0);
}