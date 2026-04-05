#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QComboBox>

struct Monom {
    int x;
    int y;
};

class PlotWidget : public QWidget {
    Q_OBJECT
public:
    explicit PlotWidget(QWidget *parent = nullptr);
    void setGenerators(const QVector<Monom>& gens);
    void setMode(bool showIdeal);
    
protected:
    void paintEvent(QPaintEvent *event) override;
    
private:
    QVector<Monom> m_generators;
    bool m_showIdeal = true;
    static constexpr int MAX_COORD = 15;
};

class MainWindow : public QMainWindow {
    Q_OBJECT
    
public:
    MainWindow(QWidget *parent = nullptr);
    
private slots:
    void onGenerateClicked();
    void onModeChanged(int index);
    
private:
    QVector<Monom> parseGenerators(const QString& text);
    
    QLineEdit* m_inputEdit;
    QPushButton* m_generateButton;
    QComboBox* m_modeCombo;
    QLabel* m_statusLabel;
    PlotWidget* m_plotWidget;
};

#endif